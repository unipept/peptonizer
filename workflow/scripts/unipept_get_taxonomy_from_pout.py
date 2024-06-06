import requests
import json
import re
import os.path
import argparse
import logging
import time
import asyncio
import aiohttp
import aiofiles
import sys
import gzip

from ete3 import NCBITaxa

UNIPEPT_URL = "http://rick.ugent.be/mpa/pept2filtered.json"

parser = argparse.ArgumentParser()

parser.add_argument(
    "--unipept-response-file",
    type=str,
    required=True,
    help="Output: path to Unipept response .json file",
)
parser.add_argument(
    "--taxonomy-query",
    required=True,
    help="Taxa that should be used to query in Unipept. If querying all taxa, put [1].",
)
parser.add_argument(
    "--fdr",
    type=float,
    required=True,
    help="Min peptide score for the peptide to be included in the search.",
)
parser.add_argument(
    "--pout-file",
    type=str,
    nargs="+",
    required=True,
    help="Input: paths to percolator (ms2rescore) '.pout' files.",
)
parser.add_argument("--unipept-peptide-counts", type=str, required=True, help="Path to output file that contains all queried peptide counts (which should be used in the next step).")
parser.add_argument(
    "--log-file",
    type=str,
    required=True,
    help="Output: path to logfile where failed Unipept query attempts are stored.",
)

args = parser.parse_args()

ncbi = NCBITaxa()


def Poutparser(pout_files, fdr_threshold, decoy_flag):
    """
    Parses the ms2rescore pout file for peptides, psm numbers and peptide scores.

    Note: this code was adapted from the Pout2Prot tool.

    :param pout_files: str, path to pout file(s)
    :param fdr_threshold: float, fdr threshold below which psms are kept
    :param decoy_flag: str, can be emtpy string, decoy flag in pout file
    :return: dict, peptides:[score,#psms]
    """

    pep_score = dict()
    pep_psm = dict()
    pep_score_psm = dict()

    for pout_file in pout_files:
        with gzip.open(pout_file, "rt") as f:
            next(f)  # skip header
            for line in f:
                # skip empty lines
                if line.rstrip() == "":
                    continue
                splitted_line = line.rstrip().split("\t", maxsplit=5)
                assert (
                    len(splitted_line) >= 6
                ), "Input file is wrongly formatted. Make sure that the input is a valid .pout file."
                psm_id, _, q, pep, peptide, _ = splitted_line
                if float(q) < fdr_threshold:
                    peptide = re.sub("\[.*?\]", "", peptide)
                    peptide = peptide.split(".")[1]
                    # update pep_psm
                    if peptide not in pep_psm.keys():
                        pep_psm[peptide] = set()
                        pep_psm[peptide].add(psm_id)
                    else:
                        pep_psm[peptide].add(psm_id)
                    # update pep_score
                    if peptide not in pep_score.keys():
                        if float(pep) < 0.001:
                            pep_score[peptide] = "0.001"
                        else:
                            pep_score[peptide] = (
                                pep  # adjustement necessary to not have 0 and 1 fuck up probability calculations
                            )
                    else:
                        if float(pep) < 0.001:
                            pep_score[peptide] = "0.001"
                        else:
                            pep_score[peptide] = min(pep, pep_score[peptide])
                    pep_score_psm[peptide] = [pep_score[peptide], len(pep_psm[peptide])]

    return pep_score_psm


def MS2RescoreOutParser(pout_files, fdr_threshold, decoy_flag):
    """
    Parses the ms2rescore pout file for peptides, psm numbers and peptide scores
    :param pout_file: str, path to pout file(s)
    :param fdr_threshold: float, fdr threshold below which psms are kept
    :param decoy_flag: str, can be emtpy string, decoy flag in pout file
    :return: dict, peptides:[score,#psms]
    """

    pep_score = dict()
    pep_psm = dict()
    pep_score_psm = dict()

    for pout_file in pout_files:
        with gzip.open(pout_file, "rt") as f:
            next(f)  # skip header
            for line in f:
                # skip empty lines
                if line.rstrip() == "":
                    continue
                splitted_line = line.rstrip().split("\t")[0:8]
                # assert len(splitted_line) >= 6, "Input file is wrongly formatted. Make sure that the input is a valid .pout file."
                peptide, psm_id, run, colelction, collection, score, q, pep = (
                    splitted_line
                )
                if float(q) < fdr_threshold:
                    peptide = re.sub("\[.*?\]", "", peptide)
                    peptide = peptide.split("/")[0]
                    # update pep_psm
                    if peptide not in pep_psm.keys():
                        pep_psm[peptide] = set()
                        pep_psm[peptide].add(psm_id)
                    else:
                        pep_psm[peptide].add(psm_id)
                    # update pep_score
                    if peptide not in pep_score.keys():
                        if float(pep) < 0.001:
                            pep_score[peptide] = "0.001"
                        else:
                            pep_score[peptide] = (
                                pep  # adjustement necessary to not have 0 and 1 fuck up probability calculations
                            )
                    else:
                        if float(pep) < 0.001:
                            pep_score[peptide] = "0.001"
                        else:
                            pep_score[peptide] = min(pep, pep_score[peptide])
                    pep_score_psm[peptide] = [pep_score[peptide], len(pep_psm[peptide])]

    return pep_score_psm


def PepListNoMissedCleavages(peptide):
    """
    Takes a peptide and cleaves it into Unipept format (0 missed cleavages,
    cleaves after K or R except followed by P)
    :param peptides_in: list of peptides
    :return: list of peptides in Unipept format
    """

    peptides = list()
    trypsin = lambda peptide: re.sub(
        r"(?<=[RK])(?=[^P])", "\n", peptide, re.DOTALL
    ).split()
    peptides += trypsin(peptide)

    return peptides


def generatePostRequestChunks(peptides, TargetTaxa, chunksize=1000, cutoff=1000):
    """
    Generates POST requests (json) querying a chunk of peptides from peptide list and target taxon
    :param peptides: list of peptides to query in Unipept
    :param TargetTaxa: list of one or more taxa to include in the Unipept query
    :param chunksize: number of peptides to be requested from Unipept
    :param cutoff: number of proteins a peptide is associated to above which said peptide will be removed from the query by Unipept. This enhances query speed.
    """
    print("querying taxa ", TargetTaxa)
    AllTargetTaxa = []
    for Taxon in TargetTaxa:
        AllTargetTaxa.append(Taxon)
        AllTargetTaxa.extend(ncbi.get_descendant_taxa(Taxon, collapse_subspecies=True))

    Listofpeptides = [
        peptides[i : i + chunksize] for i in range(0, len(peptides), chunksize)
    ]
    Listofrequests = [
        {"cutoff": cutoff, "peptides": chunk, "taxa": AllTargetTaxa}
        for chunk in Listofpeptides
    ]

    return Listofrequests


def PostInfoFromUnipeptChunks(request_json, out_file, failed_requests_file):
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


async def main(request, out_file, failed_requests_file):
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


pep_score_psm = MS2RescoreOutParser(args.pout_file, args.fdr, "")
UnipeptPeptides = dict()
for peptide in pep_score_psm.keys():
    FullyTrypticPeptides = PepListNoMissedCleavages(peptide)
    for pep in FullyTrypticPeptides:
        UnipeptPeptides[pep] = {
            "score": pep_score_psm[peptide][0],
            "psms": pep_score_psm[peptide][1],
        }
with open(args.unipept_peptide_counts, "a") as f_out:
    f_out.write(json.dumps(UnipeptPeptides))

# get and save Info from Unipept if the response file doesn't exist yet
request = generatePostRequestChunks(
    list(UnipeptPeptides.keys()), [int(item) for item in args.taxonomy_query.split(",")]
)
save = asyncio.run(main(request, args.unipept_response_file, args.log_file))
