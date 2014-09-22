import glob
import os
import re
import jacquard_utils

def find_somatic_positions(in_files, output_dir):
    somatic_positions = {}
    no_jq_tags = []
    
    total_number_of_files = len(in_files)
    count = 1
    for file in in_files:
        print "Reading [{0}] ({1}/{2})".format(os.path.basename(file), count, total_number_of_files)
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
                    format_sample_dict = jacquard_utils.combine_format_values(format_col, sample_col)
                    for key in format_sample_dict.keys():
                        if re.search(jacquard_utils.jq_filter_tag, key):
                            if format_sample_dict[key] == "1":
                                somatic_key = "^".join([split_line[0], split_line[1]])
                                somatic_positions[somatic_key] = 1
                                somatic = 1

        if somatic == 0:
            no_jq_tags.append(file)
            print "ERROR: input file [{0}] has no Jaquard tags.".format(os.path.basename(file))
            
        in_file.close()
        
        count += 1
        
    if no_jq_tags != []:
        print "ERROR: [{0}/{1}] VCF files have no Jacquard tags. Review input and try again.".format(len(no_jq_tags), len(in_files))
        exit(1)
        
    print "Found [{0}] high-confidence somatic positions".format(len(somatic_positions.keys()))
    somatic_positions_header = "##jacquard.filterHCSomatic.total_highConfidence_somatic_positions={0}\n".format(len(somatic_positions.keys()))
    print somatic_positions_header
    
    return somatic_positions, somatic_positions_header

def write_somatic(in_files, output_dir, somatic_positions, execution_context):
    for file in in_files:
        non_somatic = 0
        headers = []
        variants = []
        
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
                    variants.append(line)
                else: 
                    non_somatic += 1
       
        excluded_variants = "##jacquard.filterHCSomatic.excluded_variants={0}\n".format(non_somatic)
        print os.path.basename(file) + ": " + excluded_variants
        headers.append(excluded_variants)
        headers.extend(execution_context)
        
        sorted_headers = jacquard_utils.sort_headers(headers)
        jacquard_utils.write_output(out_file, sorted_headers, variants)
        
        in_file.close()
        out_file.close()
        
    print "Jacquard wrote [{0}] VCF files to [{1}]".format(len(in_files), output_dir)
        
    return excluded_variants
    
def filter_somatic_positions(input_dir, output_dir, execution_context=[]):
    in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
    if len(in_files) < 1:
        print "Error: Specified input directory [{0}] contains no VCF files. Check parameters and try again."
        exit(1)
        
    print "\n".join(execution_context)
    print "Processing [{0}] VCF file(s) from [{1}]".format(len(in_files), input_dir)
    
    
    somatic_positions, somatic_positions_header = find_somatic_positions(in_files, output_dir)
    execution_context.append(somatic_positions_header)
    
    write_somatic(in_files, output_dir, somatic_positions, execution_context)

def add_subparser(subparser):
    parser_normalize_vs = subparser.add_parser("filter_hc_somatic", help="Accepts a directory of Jacquard-tagged VCF results from one or more callers and creates a new directory of VCFs, where rows have been filtered to contain only positions that were called high-confidence somatic in any VCF.")
    parser_normalize_vs.add_argument("input_dir", help="Path to directory containing VCFs. All VCFs in this directory must have Jacquard-specific tags (see jacquard.py tag for more info")
    parser_normalize_vs.add_argument("output_dir", help="Path to output directory. Will create if doesn't exist and will overwrite files in output directory as necessary")

def execute(args, execution_context): 
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    
    jacquard_utils.validate_directories(input_dir, output_dir)
    
    filter_somatic_positions(input_dir, output_dir, execution_context)
    