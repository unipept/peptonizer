import requests
import json
import re
import logging
import time
import asyncio
import aiohttp
import aiofiles

from ete3 import NCBITaxa
from typing import Dict, Tuple


UNIPEPT_URL = "http://rick.ugent.be/mpa/pept2filtered.json"


def cleave_peptides(peptide: str):
    """
    Takes a peptide and cleaves it into Unipept format (0 missed cleavages, cleaves after K or R except followed by P)
    :param peptide: list of peptides
    :return: list of peptides in Unipept format
    """

    peptides = list()
    trypsin = lambda peptide: re.sub(
        r"(?<=[RK])(?=[^P])", "\n", peptide, re.DOTALL
    ).split()
    peptides += trypsin(peptide)

    return peptides


def generate_post_request_chunks(peptides, target_taxa, chunk_size=1000, cutoff=1000):
    """
    Generates POST requests (json) querying a chunk of peptides from peptide list and target taxon
    :param peptides: list of peptides to query in Unipept
    :param target_taxa: list of one or more taxa to include in the Unipept query
    :param chunk_size: number of peptides to be requested from Unipept
    :param cutoff: number of proteins a peptide is associated to above which said peptide will be removed from the query by Unipept. This enhances query speed.
    """
    print("querying taxa ", target_taxa)

    ncbi = NCBITaxa()

    all_target_taxa = []
    for Taxon in target_taxa:
        all_target_taxa.append(Taxon)
        all_target_taxa.extend(ncbi.get_descendant_taxa(Taxon, collapse_subspecies=True))

    list_of_peptides = [
        peptides[i: i + chunk_size] for i in range(0, len(peptides), chunk_size)
    ]

    list_of_requests = [
        {"cutoff": cutoff, "peptides": chunk, "taxa": all_target_taxa}
        for chunk in list_of_peptides
    ]

    return list_of_requests


def post_info_from_unipept_chunks(request_json, out_file, failed_requests_file):
    """
    Send all requests, get for each peptide the phylum, family, genus and collection of EC-numbers
    :param request_list: list of Get Requests
    :param result_file: csv file with Unipept info (phylum, family, genus and collection of EC-numbers)
    :return: None
    """
    logging.basicConfig(filename=failed_requests_file, level=logging.INFO)

    print("now querying Unipept in " + str(len(request_json)) + " chunks")

    # try Unipept query in chunks
    failed_requests = {}
    count = 0
    for chunk in request_json:
        start = time.time()
        try:

            request = requests.post(
                UNIPEPT_URL,
                json.dumps(chunk),
                headers={"content-type": "application/json"},
                timeout=1800,
            )
            request.raise_for_status()

            with open(out_file, "a") as f_out:
                print(request.text, file=f_out)
            end = time.time()
            print(
                "sucessfully queried a chunk, number "
                + str(count)
                + " out of "
                + str(len(request_json))
                + "\n time taken"
                + str(end - start)
            )
            count += 1

        except requests.exceptions.RequestException as e:

            logging.error(f"Request {UNIPEPT_URL} failed with error {e}")
            failed_requests[json.dumps(chunk)] = e

    # retry failed requests
    for chunk, error in failed_requests.items():
        try:
            request = requests.post(
                UNIPEPT_URL, chunk, headers={"content-type": "application/json"}, timeout=3600
            )
            request.raise_for_status()
            with open(out_file, "a") as f_out:
                print(request.text, file=f_out)

        except requests.exceptions.RequestException as e:
            logging.error(f"Retry request to {UNIPEPT_URL} failed with error: {e}")


async def fetch_data(
    session, url, json_input, out_file, failed_requests, i, total_chunks
):
    try:
        async with session.post(
            url,
            json=json_input,
            timeout=8600,
            headers={"content-type": "application/json"},
        ) as response:
            response.raise_for_status()
            result = await response.text()
            async with aiofiles.open(out_file, "a") as f_out:
                await f_out.write(result + "\n")
            print(f"Successfully queried chunk {i} out of {total_chunks}")
            return True
    except (
        aiohttp.ClientError,
        aiohttp.ClientConnectionError,
        asyncio.TimeoutError,
    ) as e:
        print(f"Failed to query chunk {i}")
        logging.error(f"Request to {url} failed with error: {e}")
        failed_requests.append((json_input, e))
        return False


async def limited_gather(
    semaphore, session, url, chunks, out_file, failed_requests, i, total_chunks
):
    async with semaphore:
        return await fetch_data(
            session, url, chunks, out_file, failed_requests, i, total_chunks
        )


async def query(request, out_file, failed_requests_file):
    logging.basicConfig(filename=failed_requests_file, level=logging.INFO)
    print(f"Now querying Unipept in {len(request)} chunks")

    semaphore = asyncio.Semaphore(3)  # Limit the number of concurrent requests to three

    total_chunks = len(request)
    failed_requests = []

    async with aiohttp.ClientSession() as session:
        limited_gather_tasks = []
        for i, chunk in enumerate(request):
            task = limited_gather(
                semaphore,
                session,
                UNIPEPT_URL,
                chunk,
                out_file,
                failed_requests,
                i,
                total_chunks,
            )
            limited_gather_tasks.append(task)

        await asyncio.gather(*limited_gather_tasks)

        # Retry failed requests with exponential backoff
        for retry in range(3):  # Retry a maximum of 3 times
            retry_tasks = []
            for chunk, e in failed_requests:
                logging.error(f"Request {chunk} try 2 to {UNIPEPT_URL} failed with error: {e}")
                task = limited_gather(
                    semaphore,
                    session,
                    UNIPEPT_URL,
                    chunk,
                    out_file,
                    failed_requests,
                    i,
                    total_chunks,
                )
                retry_tasks.append(task)

            await asyncio.gather(*retry_tasks)
            await asyncio.sleep(2**retry)  # Exponential backoff


def fetch_unipept_taxon_information(
    pep_score_psm: Dict[str, Tuple[float, int]],
    unipept_peptide_counts_file: str,
    unipept_response_file: str,
    taxonomy_query: str,
    log_file: str
):
    unipept_peptides = dict()

    for peptide in pep_score_psm.keys():
        fully_tryptic_peptides = cleave_peptides(peptide)
        for pep in fully_tryptic_peptides:
            unipept_peptides[pep] = {
                "score": pep_score_psm[peptide][0],
                "psms": pep_score_psm[peptide][1],
            }

    with open(unipept_peptide_counts_file, "a") as f_out:
        f_out.write(json.dumps(unipept_peptides))

    # Get and save info from Unipept if the response file doesn't exist yet
    request = generate_post_request_chunks(
        list(unipept_peptides.keys()), [int(item) for item in taxonomy_query.split(",")]
    )

    asyncio.run(query(request, unipept_response_file, log_file))
