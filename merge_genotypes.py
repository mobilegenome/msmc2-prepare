"""
Author: Fritjof Lammers
Date: 2023-11-14

License: MIT

This script reads in a VCF file for a single sample (called "A") and another VCF-like file that contains phased genotypes (called "B"). 

The unphased genotype in file A is indicated by "0/0", "0/1", or "1/1". 
The phased genotype in file B is indicated by "0|0", "0|1", or "1|1". 

The script then iterates over genomic coordinates in both files and does the following operations: 
 1. If the coordinate is present in both files, the genotype from file B is used to replaced the unphased genotype in file A. 
 2. If the coordinate is present in file A but not in file B, the genotype from file A is used, i.e. left unchanged.
 3. If the coordinate is present in file B but not in file A, the genotype is skipped. 

In all cases where matching coordinates in file A and B are present the REF and ALT allele are checked for consistency.

All remaining contents from file A are kept to ensure adherence to the VCF format.
"""

import sys
import gzip

def phase_genotype(gt): 
    return gt.replace('/', '|')

ALLOWED_GENOTYPES = {'0/0', '1/0', '0/1', '1/1'}

EVENT_LOGGER = {
        "read_lines":                      [],
        "inconsistent_ref_alt_genotypes":  [],
        "inconsistent_ref_alt_alleles":    [],
        "replaced_genotypes":              [],
        "kept_genotypes":                  [],
    }

VERBOSE = False

def open_file(file_path, _mode='rt'):
    """
    Opens a file, handling both regular and gzipped files.

    Args:
    - file_path (str): The path to the file.

    Returns:
    - file object: An open file object.
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, _mode)
    else:
        return open(file_path, 'r')

def read_vcf_header(file_path):
    """
    Reads a VCF file and returns list with the header lines

    Args:
    - file_path (str): The path to the VCF file.

    Returns:
    - list: A list with header lines
    """
    vcf_data = []
    with open_file(file_path) as file:
        for line in file:
            if line.startswith('#'):
                vcf_data.append(line)
    return vcf_data


def read_vcf(file_path):
    """
    Reads a VCF file and returns a dictionary where the keys are genomic coordinates
    and the values are the corresponding genotypes.

    Args:
    - file_path (str): The path to the VCF file.

    Returns:
    - dict: A dictionary with genomic coordinates as keys and genotypes as values.
    """
    vcf_data = {}
    with open_file(file_path) as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                coordinate = (parts[0], int(parts[1]))
                ref = parts[3]
                alt = parts[4]
                format_genotype = parts[9].split(':')[0]
                format_other = "".join(parts[9].split(':')[1:])
                vcf_data[coordinate] = (line, ref, alt, format_genotype, format_other)
    return vcf_data

def replace_genotypes(file_a, file_b):
    """
    Replaces unphased genotypes in file A with phased genotypes from file B based on genomic coordinates.

    Args:
    - file_a (str): Path to the VCF file for sample A.
    - file_b (str): Path to the VCF-like file containing phased genotypes.

    Returns:
    - str: The updated content for sample A in VCF format.
    """
    vcf_a = read_vcf(file_a)
    vcf_b = read_vcf(file_b)

    updated_content = []

    for coordinate, (line_a, ref_a, alt_a, format_genotype_a, format_other_a) in vcf_a.items():
        EVENT_LOGGER["read_lines"].append(coordinate)
        if coordinate in vcf_b:
            (line_b, ref_b, alt_b, format_genotype_b, format_other_b) = vcf_b[coordinate]

            # Check consistency of REF and ALT alleles
            if format_genotype_a not in ALLOWED_GENOTYPES:
                sys.exit("Error: Unphased genotype in file A must be {",".join(ALLOWED_GENOTYPES)}")
            if format_genotype_b not in map(phase_genotype, ALLOWED_GENOTYPES):
                sys.exit("Error: Phased genotype in file B must be {",".join(map(phase_genotype, ALLOWED_GENOTYPES))}")

            ref_alt_a = format_genotype_a.split('/')
            ref_alt_b = format_genotype_b.split('|')

            if ref_alt_a != ref_alt_b:
                EVENT_LOGGER["inconsistent_ref_alt_genotypes"].append(coordinate)

                if VERBOSE:
                    print(f"INFO: Inconsistent REF and ALT genotypes at coordinate {coordinate}: {ref_alt_a} vs {ref_alt_b}", file=sys.stderr)
            
            alleles_a = ref_a, alt_a
            alleles_b = ref_b, alt_b

            if alleles_a != alleles_b:
                EVENT_LOGGER["inconsistent_ref_alt_alleles"].append(coordinate)

                if VERBOSE:
                    print(f"WARNING: Inconsistent REF and ALT alleles at coordinate {coordinate}: {alleles_a} vs {alleles_b}", file=sys.stderr)

            
            # Inejct genotype from file B into line from file A
            line_b_parts = line_b.strip().split('\t')
            modified_line_a_parts = line_a.strip().split('\t')
            modified_line_a_parts[9] = f"{format_genotype_b}:{format_other_a}"
            modified_line_a = '\t'.join(modified_line_a_parts) + '\n'
            updated_content.append(modified_line_a)
            EVENT_LOGGER["replaced_genotypes"].append(coordinate)

            if VERBOSE:
                print("INFO: Replaced genotype at coordinate", coordinate, file=sys.stderr)
        else:
            # Coordinate not present in file B, use genotype from file A
            updated_content.append(line_a)
            EVENT_LOGGER["kept_genotypes"].append(coordinate)

            if VERBOSE:
                print("INFO: Kept genotype at coordinate", coordinate, file=sys.stderr)
            


    # Add remaining contents from file A
    # for coordinate, (line_a,ref_a, alt_a, format_gentype_a, format_other_a, _) in vcf_a.items():
    #     if coordinate not in vcf_b:
    #         updated_content.append(line_a)
    #         print("INFO: Kept genotype at coordinate", coordinate, file=sys.stderr)


    return ''.join(updated_content)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: python script.py <file_A.vcf[.gz]> <file_B.vcf[.gz]> <output.vcf[.gz]")
    
    file_a_path = sys.argv[1]
    file_b_path = sys.argv[2]
    output_file_path = sys.argv[3]

    header_a  = read_vcf_header(file_a_path)
    result = replace_genotypes(file_a_path, file_b_path)
    with open_file(output_file_path, "wt") as output_file:
        output_file.write("".join(header_a))
        output_file.write(result)
    
    # print event counts in EVENT_LOGGER
    for event, coordinates in EVENT_LOGGER.items():
        print(f"{event}: {len(coordinates)}")

