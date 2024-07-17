#!/bin/env python3

import argparse
import gzip
import os
from glob import glob

class LINE:
    def __init__(self):
        self.phenotype  = {}

    def fill(self, phenotypes: dict):
        for pheno, is_float in phenotypes.items():
            if (pheno not in self.phenotype):
                if (is_float):
                    self.phenotype[pheno] = -9.0
                else:
                    self.phenotype[pheno] = -9
    
    def get_phenotypes(self):
        phenos = []
        for p in self.phenotype.keys():
            phenos.append(p)
        phenos.sort()
        return phenos
  
### functions ###
def get_arguments():
    """get the input files"""

    parser = argparse.ArgumentParser(
        description="Reformat input files into ped and pheno maps")
    parser.add_argument("-p", "--pheno_dir", help="Path to phenotype csv files", required=True)
    parser.add_argument("-v", "--vcf", help="Input vcf file containing snp information", required=True)
    parser.add_argument("--no_ped", help="If set, do not write ped file", action="store_true")

    args       = parser.parse_args()
    pheno_dir  = args.pheno_dir
    vcf        = args.vcf
    no_ped     = args.no_ped

    assert os.path.isfile(vcf), f"Could not locate {vcf}"
    assert os.path.isdir(pheno_dir), f"Could not locate {pheno_dir}"

    return vcf, pheno_dir, no_ped

def get_genotypes(ref: str, alt: str, genotype: str):
    """return the genotype of the sample"""

    if (genotype[0] == '.'): return (0,0)

    c = (ref, alt)

    a1, a2 = genotype.split('/')
    a1, a2 = int(a1), int(a2)
    a1, a2 = c[a1], c[a2]

    return (a1, a2)


def write_ped(vcf: str, outfile: str, no_ped):
    """parse through vcf to create the ped file"""

    snp_map = {}
    snps    = []

    # open up the file stream
    fh = gzip.open(vcf, 'rt') if vcf.endswith(".gz") else open(vcf, 'r')

    for line in fh:
        # skip empty
        if (len(line) == 0): continue
        if (line[0] == '#'):
            if (line[1] == '#'): continue
            fields  = line.strip('\n').split('\t')
            samples = {}
            samples_list = []
            for i in range(9, len(fields)):
                samples[i]         = fields[i]
                samples_list.append(fields[i])
                snp_map[fields[i]] = []
            if (no_ped):
                break
            else:
                continue
        fields = line.strip('\n').split('\t')
        pos    = int(fields[1])
        chr    = fields[0]
        ref    = fields[3]
        alt    = fields[4]
        snps.append(f"{chr}_{pos}")
        for i in range(9, len(fields)):
            genotype = get_genotypes(ref, alt, fields[i])
            snp_map[samples[i]].append(genotype)
    fh.close()

    if (no_ped):
        print(f"Found {len(samples_list)} samples")
        return samples_list

    # write .ped file
    cnt    = 0
    ofh    = open(outfile, 'w')
    header = "FID\tID\tPID\tMID\tSex\tPhenotype\t" + '\t'.join(snps) + '\n'
    ofh.write(header)
    for sample, genos in snp_map.items():
        cnt   += 1
        start  = f"{cnt}\t{sample}\t0\t0\t0\t-9"
        end    = ''
        for geno in genos:
            end += f"\t{geno[0]} {geno[1]}" 
        outline = start + end + '\n'
        ofh.write(outline)

    ofh.close()

    print("Wrote", outfile)

    # write .map file
    ofh = open(outfile.replace(".ped", ".map"), 'w')

    for snp in snps:
        fields  = snp.split('_')
        chr     = ''.join(fields[:-1])
        pos     = int(fields[-1])
        outline = f"{chr}\t{snp}\t{pos}\n"
        ofh.write(outline)
    ofh.close()

    print("Wrote", outfile.replace(".ped", ".map"))


    print(f"Found {len(samples_list)} samples")
    return samples_list

def make_outname(vcf: str, new_ext: str):
    """create the name of the output file"""

    bname = os.path.basename(vcf)

    ext   = ".vcf.gz" if bname.endswith(".gz") else ".vcf"
    oname = bname.replace(ext, new_ext)

    return oname

def get_pheno(input_file: str):
    """get phenotype trait based on input file"""

    bname      = os.path.basename(input_file)
    ext        = ".female.csv" if bname.endswith(".female.csv") else ".csv"
    pheno_type = bname.replace(ext, '')

    return pheno_type

def parse_phenos(pheno_dir: str, samples: list):
    """parse each phenotype file for the values"""

    samples_map = {s : LINE() for s in samples} # Sample : {Phenotype: N}
    pheno_map   = {}
    if (pheno_dir[-1] != '/'):
        pheno_dir += '/'

    pheno_files = glob(f"{pheno_dir}*.csv")

    for pfile in pheno_files:
        phenotype = get_pheno(pfile)
        fh = gzip.open(pfile, "rt") if pfile.endswith(".gz") else open(pfile, 'r')
        print("Parsing", pfile)
        for line in fh:
            fields = line.strip('\n').split(',')
            sample = fields[0]
            if ('.' in fields[1]):
                pheno_val = float(fields[1])
                if (phenotype not in pheno_map):
                    pheno_map[phenotype] = True
            else:
                try:
                    pheno_val = int(fields[1])
                except:
                    print("Found a weird line:", line)
                    continue
            if (sample not in samples_map):
                print(f"{sample} not in vcf file")
                continue
            samples_map[sample].phenotype[phenotype]= pheno_val
        if (phenotype not in pheno_map):
            pheno_map[phenotype] = False
        fh.close()
        print("Finished", pfile)

    for line_obj in samples_map.values():
        line_obj.fill(pheno_map)

    return samples_map

def write_phenos(samples: list, samples_map: dict, outname: str):
    """write the phenotype file using the map of phenotypes"""

    # create the phenotype output
    ofh        = open(outname, 'w')
    phenotypes = None

    for i, sample in enumerate(samples, start=1):
        line_obj = samples_map[sample]
        if (i == 1):
            phenotypes = line_obj.get_phenotypes()
            header = f"FID\tID\t" + '\t'.join(phenotypes) + '\n'
            ofh.write(header)
        values  = [str(line_obj.phenotype[p]) for p in phenotypes]
        outline = f"{i}\t{sample}\t" + '\t'.join(values) + '\n'
        ofh.write(outline)

    print("Wrote", outname)

    ofh.close()

def main():
    """start the pipeline"""

    #get files
    vcf, pheno_dir, no_ped = get_arguments()

    # create the output name
    outfile = make_outname(vcf, ".ped")
    
    # parse vcf & write ped
    samples = write_ped(vcf, outfile, no_ped)

    # parse csv files and create .pheno file
    samples_map = parse_phenos(pheno_dir, samples)

    outfile = make_outname(vcf, ".pheno")

    # write out phenotype file
    write_phenos(samples, samples_map, outfile)

if __name__ == "__main__":
    main()

