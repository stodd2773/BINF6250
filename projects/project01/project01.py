#!/usr/bin/env python
"""
Authors: Spencer Todd, Meghana Ravi, Nicholas Bottomley

Reads a clinvar vcf file and identifies diseases
associated with rare allele variants.

Prints a list of the frequency of rare variant
associated diseases.

Usage: python3 project01.py
"""


from pprint import pprint


def parse_line(line):
    """
    Parses a VCF line to extract rare variant diseases.

    Args:
        line (str): A single line from a VCF file

    Returns:
        list: A list of diseases associated with rare variants
    """
    columns = line.split('\t')
    if len(columns) < 8:
        return []
    info_column = columns[7]
    info_parts = info_column.split(';')

    info_dict = {}
    for part in info_parts:
        if '=' in part:
            key, value = part.split('=')
            info_dict[key] = value

    af_exac = info_dict.get('AF_EXAC')
    if af_exac is None:
        return []

    af_exac_value = float(af_exac)

    if af_exac_value < 0.0001:
        clndn = info_dict.get('CLNDN')
        diseases = [disease for disease in clndn.split('|')
                    if disease not in ('not_specified', 'not_provided')]
    else:
        return []

    return diseases


def read_file(file_name):
    """
    Takes in file name string argument,
    then opens and reads line by line,
    passing line to parsing function to
    increment count of rare-variant associated
    diseases in given file.

    Args:
        file_name (str): filepath for vcf file to be parsed

    Returns:
        dict: Dictionary containing count of disease frequencies
    """
    disease_count = {}

    with open(file_name, mode='r', encoding='utf-8') as infile:
        for line in infile:
            clean_line = line.rstrip()

            # Skip information and header lines
            if clean_line[0] == "#":
                continue

            disease_parse = parse_line(line)

            # If parse_line returns one or more disease
            if len(disease_parse) >= 1:
                for disease in disease_parse:
                    if disease in disease_count:
                        disease_count[disease] += 1
                    else:
                        disease_count[disease] = 1

    return disease_count


if __name__ == "__main__":
    pprint(read_file("clinvar_20190923_short.vcf"))
