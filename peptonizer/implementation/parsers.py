import gzip
import re

from typing import Dict, Tuple, List


def parse_pout(pout_files: List[str], fdr_threshold: float, decoy_flag: str) -> Dict[str, Tuple[float, int]]:
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
                                pep  # adjustment necessary to not have 0 and 1 fuck up probability calculations
                            )
                    else:
                        if float(pep) < 0.001:
                            pep_score[peptide] = "0.001"
                        else:
                            pep_score[peptide] = min(pep, pep_score[peptide])
                    pep_score_psm[peptide] = [pep_score[peptide], len(pep_psm[peptide])]

    return pep_score_psm


def parse_ms2rescore_output(pout_files: List[str], fdr_threshold: float, decoy_flag: str) -> Dict[str, Tuple[float, int]]:
    """
    Parses the ms2rescore pout file for peptides, psm numbers and peptide scores
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
                splitted_line = line.rstrip().split("\t")[0:8]
                # assert len(splitted_line) >= 6, "Input file is wrongly formatted. Make sure that the input is a valid .pout file."
                peptide, psm_id, run, collection, is_decoy, score, q, pep = (
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
                                pep  # adjustment necessary to not have 0 and 1 fuck up probability calculations
                            )
                    else:
                        if float(pep) < 0.001:
                            pep_score[peptide] = "0.001"
                        else:
                            pep_score[peptide] = min(pep, pep_score[peptide])
                    pep_score_psm[peptide] = [pep_score[peptide], len(pep_psm[peptide])]

    return pep_score_psm