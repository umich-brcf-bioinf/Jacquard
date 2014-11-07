import glob
import os
import re
import utils

import logger as logger

def find_somatic_positions(in_files, output_dir):
    somatic_positions = {}
    no_jq_tags = []
    
    total_number_of_files = len(in_files)
    count = 1
    for file in in_files:
        logger.info("Reading [{}] ({}/{})", os.path.basename(file), count, total_number_of_files)
        somatic = 0
        in_file = open(file, "r")
        
        for line in in_file:
            if line.startswith("#"):
                continue
            else:
                split_line = line.split("\t")
                format_col = split_line[8]
                sample_cols = split_line[9:]
                
                for sample_col in sample_cols:
                    format_sample_dict = utils.combine_format_values(format_col, sample_col)
                    for key in format_sample_dict.keys():
                        if re.search(utils.jq_somatic_tag, key):
                            if format_sample_dict[key] == "1":
                                somatic_key = "^".join([split_line[0], split_line[1]])
                                somatic_positions[somatic_key] = 1
                                somatic = 1

        if somatic == 0:
            no_jq_tags.append(file)
            logger.warning("Input file [{}] has no high-confidence somatic variants.", os.path.basename(file))
            
        in_file.close()
        
        count += 1
        
    if no_jq_tags:
        logger.error("[{}/{}] VCF files had no high-confidence somatic variants. Review input and try again.", len(no_jq_tags), len(in_files))
        exit(1)
        
    logger.info("Found [{}] high-confidence somatic positions", len(somatic_positions.keys()))
    somatic_positions_header = "##jacquard.filterHCSomatic.total_highConfidence_somatic_positions={0}\n".format(len(somatic_positions.keys()))
#     logger.info("{}", somatic_positions_header)
    
    return somatic_positions, somatic_positions_header

def write_somatic(in_files, output_dir, somatic_positions, execution_context):
    for file in in_files:
        non_somatic = 0
        headers = []
        actual_sorted_variants = []
        
        fname, extension = os.path.splitext(os.path.basename(file))
        new_file = fname + "_HCsomatic" + extension
        in_file = open(file, "r")
        
        out_file = open(os.path.join(output_dir, new_file), "w")
        
        for line in in_file:
            if line.startswith("#"):
                headers.append(line)
            else:
                split_line = line.split("\t")
                key = "^".join([split_line[0], split_line[1]])
                if key in somatic_positions:
                    actual_sorted_variants.append(line)
                else: 
                    non_somatic += 1
       
        excluded_variants = "##jacquard.filterHCSomatic.excluded_variants={0}\n".format(non_somatic)
#         logger.info("{}:{}", os.path.basename(file), excluded_variants)
        headers.append(excluded_variants)
        headers.extend(execution_context)
        
        sorted_headers = utils.sort_headers(headers)
        utils.write_output(out_file, sorted_headers, actual_sorted_variants)
        
        in_file.close()
        out_file.close()
        
    logger.info("Jacquard wrote [{}] VCF files to [{}]", len(in_files), output_dir)
        
    return excluded_variants
    
def filter_somatic_positions(input_dir, output_dir, execution_context=[]):
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    if len(in_files) < 1:
        logger.error("Specified input directory [{}] contains no VCF files. Check parameters and try again.", input_dir)
        exit(1)

    logger.info("Processing [{}] VCF file(s) from [{}]", len(in_files), input_dir)


    somatic_positions, somatic_positions_header = find_somatic_positions(in_files, output_dir)
    execution_context.append(somatic_positions_header)

    write_somatic(in_files, output_dir, somatic_positions, execution_context)

def add_subparser(subparser):
    parser = subparser.add_parser("filter_hc_somatic", help="Accepts a directory of Jacquard-tagged VCF results from one or more callers and creates a new directory of VCFs, where rows have been filtered to contain only positions that were called high-confidence somatic in any VCF.")
    parser.add_argument("input", help="Path to directory containing VCFs. All VCFs in this directory must have Jacquard-specific tags (see jacquard.py tag for more info")
    parser.add_argument("output", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)

    utils.validate_directories(input_dir, output_dir)

    filter_somatic_positions(input_dir, output_dir, execution_context)
