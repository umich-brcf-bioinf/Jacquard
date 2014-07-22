#!/usr/bin/python2.7
import jacquard_utils

def add_subparser(subparser):
    parser_normalize_vs = subparser.add_parser("filter_somatic", help="Accepts a directory of Jacquard-tagged VCF results from one or more callers and creates a new directory of VCFs, where rows have been filtered to contain only positions that were called somatic in any VCF.")
    parser_normalize_vs.add_argument("input_dir", help="Path to directory containing VCFs. Each sample must have exactly two files matching these patterns: <sample>.indel.vcf, <sample>.snp.vcf, <sample>.indel.Germline.hc, <sample>.snp.Germline.hc, <sample>.indel.LOH.hc, <sample>.snp.LOH.hc, <sample>.indel.Somatic.hc, <sample>.snp.somatic.hc, <sample>.indel.*.hc, <sample>.snp.*.hc ")
    parser_normalize_vs.add_argument("output_dir", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    jacquard_utils.validate_directories(input_dir, output_dir)