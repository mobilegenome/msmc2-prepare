"""
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
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                coordinate = (parts[0], int(parts[1]))
                genotype = parts[9].split(':')[0]
                vcf_data[coordinate] = genotype
    return vcf_data

def replace_genotypes(file_a, file_b):
    """
    Replaces unphased genotypes in file A with phased genotypes from file B based on genomic coordinates.

    Args:
    - file_a (str): Path to the VCF file for sample A.
    - file_b (str): Path to the VCF-like file containing phased genotypes.

    Returns:
    - dict: A dictionary with updated genotypes for sample A.
    """
    vcf_a = read_vcf(file_a)
    vcf_b = read_vcf(file_b)

    updated_genotypes = {}

    for coordinate, genotype_a in vcf_a.items():
        if coordinate in vcf_b:
            genotype_b = vcf_b[coordinate]

            # Check consistency of REF and ALT alleles
            if genotype_a != '0/0' and genotype_a != '0/1' and genotype_a != '1/1':
                sys.exit("Error: Unphased genotype in file A must be '0/0', '0/1', or '1/1'")
            if genotype_b != '0|0' and genotype_b != '0|1' and genotype_b != '1|1':
                sys.exit("Error: Phased genotype in file B must be '0|0', '0|1', or '1|1'")

            ref_alt_a = genotype_a.split('/')
            ref_alt_b = genotype_b.split('|')

            if ref_alt_a != ref_alt_b:
                sys.exit(f"Error: Inconsistent REF and ALT alleles at coordinate {coordinate}")

            updated_genotypes[coordinate] = genotype_b
        else:
            # Coordinate not present in file B, use genotype from file A
            updated_genotypes[coordinate] = genotype_a

    return updated_genotypes

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python script.py <file_A.vcf> <file_B.vcf>")
    
    file_a_path = sys.argv[1]
    file_b_path = sys.argv[2]

    result = replace_genotypes(file_a_path, file_b_path)

    # Print the updated genotypes
    for coordinate, genotype in result.items():
        print(f"{coordinate[0]}\t{coordinate[1]}\t{genotype}")
