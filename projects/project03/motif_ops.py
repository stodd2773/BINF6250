import numpy as np
from seq_ops import reverse_complement
from typing import List, Tuple, Optional


TRANS_TAB = str.maketrans('ACGT', '0123')


def build_pfm(sequences: List[str], length: int) -> np.ndarray:
    """
    Build a Position Frequency Matrix (PFM) from a list of sequences.

    Args:
        sequences (List[str]): List of sequence strings.
        length (int): Size of the PFM to build.

    Returns:
        np.ndarray: PFM with dimensions 4 x length.
    """
    pfm = np.zeros((4, length), dtype=np.int32)
    base_to_index = np.zeros(256, dtype=np.int8)
    base_to_index[ord('A')] = 0
    base_to_index[ord('C')] = 1
    base_to_index[ord('G')] = 2
    base_to_index[ord('T')] = 3

    for seq in sequences:
        seq_array = np.frombuffer(seq[:length].encode(), dtype=np.int8)
        indices = base_to_index[seq_array]
        np.add.at(pfm, (indices, np.arange(length)), 1)

    return pfm


def build_pwm(pfm: np.ndarray) -> np.ndarray:
    """
    Build a Position Weight Matrix (PWM) from a Position Frequency Matrix (PFM).

    Args:
        pfm (np.ndarray): PFM with dimensions 4 x length.

    Returns:
        np.ndarray: PWM with dimensions 4 x length.
    """
    p = 0.25
    bg = 0.25
    
    sums = np.sum(pfm, axis=0) + 4 * p
    pwm = np.log2((pfm + p) / sums[:, np.newaxis].T) - np.log2(bg)
    
    return pwm


def score_kmer(seq: str, pwm: np.ndarray) -> float:
    """
    Score a k-mer using a Position Weight Matrix (PWM).

    The k-mer length is expected to be the same as the PWM length.

    Args:
        seq (str): k-mer to score.
        pwm (np.ndarray): PWM for scoring.

    Returns:
        float: PWM score for the k-mer.

    Raises:
        ValueError: If the k-mer and PWM are different lengths.
    """
    if len(seq) != len(pwm[0]):
        raise ValueError('K-mer and PWM are different lengths!')

    indices = np.fromiter((int(nt) for nt in seq.translate(TRANS_TAB)), dtype = int)
    return np.sum(pwm[indices, np.arange(len(indices))])


def pfm_ic(pfm: np.ndarray) -> float:
    """
    Calculate the information content of a Position Frequency Matrix (PFM).

    Args:
        pfm (np.ndarray): Position Frequency Matrix with dimensions 4 x length.

    Returns:
        float: The information content of the PFM.
    """
    p = 0.25
    sums = np.sum(pfm, axis=0) + 4 * p
    fij = (pfm + p) / sums
    ic = np.sum(2 + np.sum(fij * np.log2(fij), axis=0, where=fij > 0))
    return ic
