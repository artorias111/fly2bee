#!/bin/env python3

import argparse
import gzip
import os

class Sample:
    def __init__(self, sample_id: str):
        self.name     = sample_id
        self.variants = []

class SNP_Manager:
    def __init__(self):
        self.snps   = {}
        self.chroms = set()

    def parse_line(self, line: str):
        fields         = line.strip('\n').split('\t')
        chr, pos, _    = fields[0].split('_')
        allele         = fields[2]
        ref, alt       = fields[10].split('/')
        self.snps[int(pos)] = (ref, alt, allele)
        self.chroms.add(chr)

    def in_pos(self, pos: int, chr: str):
        return (pos in self.snps) and (chr in self.chroms)
    
    def categorize(self, fields: list):
        
        categories = []
        var        = self.snps[int(fields[1])]
        print("Processing site", fields[0], fields[1])

        for i in range(9, len(fields)):
            genotype = fields[i]
            if (genotype == "./."):
                categories.append("missing")
                continue
            a1, a2   = genotype.split('/')
            # lower    = (a1 != a2)
            ref, alt = fields[3], fields[4]
            combo    = (ref, alt)
            a1, a2   = int(a1), int(a2)
            g1, g2   = combo[a1], combo[a2]
            if (g1 == alt) or (g2 == alt):
                var_cat = (f"{fields[1]} {var[a1]}|{var[a2]}")
            else:
                var_cat = ''
            categories.append(var_cat)

        del self.snps[int(fields[1])]

        return categories
    
    def empty(self):
        return (len(self.snps) == 0)
            
### functions ###
def get_arguments():
    """get the input files"""

    parser = argparse.ArgumentParser(
        description="Identify the genotypic & amino profile for fly lines")
    parser.add_argument("-v", "--vcf", help="Input vcf file", required=True)
    parser.add_argument("-s", "--snps", help="Input cvs file containing snp information", required=True)

    args = parser.parse_args()
    vcf  = args.vcf
    snps = args.snps

    for f in [vcf, snps]:
        assert os.path.isfile(f), f"Could not locate {f}"

    return vcf, snps

def parse_vcf(vcf: str, outfile: str, snp_manager: SNP_Manager):
    """parse through vcf"""

    # open up the file stream
    fh = gzip.open(vcf, 'rt') if vcf.endswith(".gz") else open(vcf, 'r')
    ofh = open(outfile, 'w')


    for line in fh:
        # skip empty
        if (len(line) == 0): continue
        if (line[0] == '#'):
            if (line[1] == '#'): continue
            fields  = line.strip('\n').split('\t')
            samples = {i: Sample(fields[i]) for i in range(9, len(fields))}
            continue
        fields = line.strip('\n').split('\t')
        pos    = int(fields[1])
        chr    = fields[0]
        if (snp_manager.in_pos(pos, chr)):
            categories = snp_manager.categorize(fields)
            for i in range(9, len(fields)):
                samples[i].variants.append(categories[i - 9])
        if (snp_manager.empty()): break

    fh.close()

    for i in range(9, len(fields)):
        sample  = samples[i]
        new_list = []
        for j in sample.variants:
            if (j != ''):
                new_list.append(j)
        outline = sample.name + '\t' + ", ".join(new_list) + '\n'
        ofh.write(outline)

    ofh.close()

def parse_var_tsv(snps: str):
    """parse a tsv for the snps of interest"""

    snp_manager = SNP_Manager()
    fh          = gzip.open(snps, "rt") if snps.endswith(".gz") else open(snps, 'r')

    for line in fh:
        if (line.startswith("Uploaded")): continue
        snp_manager.parse_line(line)

    fh.close()

    return snp_manager

def make_outname(vcf: str):
    """create the name of the output file"""

    bname = os.path.basename(vcf)

    ext   = ".vcf.gz" if bname.endswith(".gz") else ".vcf"
    oname = bname.replace(ext, "_variants.tsv")

    return oname

def main():
    """start the pipeline"""

    #get files
    vcf, snps = get_arguments()

    # get the variants
    snp_manager = parse_var_tsv(snps)

    # create the output name
    outfile = make_outname(vcf)
    
    # parse vcf & write output
    parse_vcf(vcf, outfile, snp_manager)


if __name__ == "__main__":
    main()

