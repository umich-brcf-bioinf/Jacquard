# pylint: disable=C0111
# pylint: disable-msg=W0403
from collections import OrderedDict
import numpy
import os

import utils as utils
import logger as logger

JQ_CONSENSUS_TAG = "JQ_CONS_"

def calculate_zscore(af_mean, af_std, dp_mean, dp_std, combined_dict):
    af_range = float(combined_dict[JQ_CONSENSUS_TAG + "AF_RANGE"])
    af_zscore = (af_range - af_mean)/af_std if af_std != 0.0 else 0.0

    rounded_af_zscore = roundTwoDigits([str(af_zscore)])

    dp_range = float(combined_dict[JQ_CONSENSUS_TAG + "DP_RANGE"])
    dp_zscore = (dp_range - dp_mean)/dp_std if dp_std != 0.0 else 0.0

    rounded_dp_zscore = roundTwoDigits([str(dp_zscore)])

    combined_dict[JQ_CONSENSUS_TAG + "AF_RANGE_ZSCORE"] = str(rounded_af_zscore)
    combined_dict[JQ_CONSENSUS_TAG + "DP_RANGE_ZSCORE"] = str(rounded_dp_zscore)

    return combined_dict

def iterate_file(reader, af_range, dp_range, function):
    meta_headers = []
    header = ""
    lines = []

    af_mean = sum(af_range)/len(af_range) if af_range != [] else ""
    af_std = numpy.std(af_range) if af_range != [] else ""
    dp_mean = sum(dp_range)/len(dp_range) if dp_range != [] else ""
    dp_std = numpy.std(dp_range) if dp_range != [] else ""

    for line in reader:
        if line.startswith("##"):
            meta_headers.append(line)
        elif line.startswith("#"):
            header = line
        else:
            new_line = process_line(line, af_range, af_mean, af_std, dp_range,
                                    dp_mean, dp_std, function)
            lines.append(new_line)

    return meta_headers, header, lines

# pylint: disable=C0301
def add_zscore(meta_headers, header, lines, writer, output_file, af_range, dp_range):
    rounded_mean_af = roundTwoDigits([str(sum(af_range)/len(af_range))])
    rounded_mean_dp = roundTwoDigits([str(sum(dp_range)/len(dp_range))])
    rounded_std_af = roundTwoDigits([str(numpy.std(af_range))])
    rounded_std_dp = roundTwoDigits([str(numpy.std(dp_range))])

    consensus_meta_headers = ['##FORMAT=<ID={0}AF_RANGE_ZSCORE,Number=A,Type=Integer,Description="Jacquard measure of consistency of allele frequencies among callers = (sample AF range - population mean AF range)/standard dev(population AF range)">\n'.format(JQ_CONSENSUS_TAG),
                              '##jacquard.consensus.{0}AF_RANGE_ZSCORE.mean_AF_range={1}\n'.format(JQ_CONSENSUS_TAG, rounded_mean_af),
                              '##jacquard.consensus.{0}AF_RANGE_ZSCORE.standard_deviation={1}\n'.format(JQ_CONSENSUS_TAG, rounded_std_af),
                              '##FORMAT=<ID={0}DP_RANGE_ZSCORE,Number=A,Type=Integer,Description="Jacquard measure of consistency of depth among callers = (sample DP range - population mean DP range)/standard dev(population DP range)">\n'.format(JQ_CONSENSUS_TAG),
                              '##jacquard.consensus.{0}DP_RANGE_ZSCORE.mean_DP_range={1}\n'.format(JQ_CONSENSUS_TAG, rounded_mean_dp),
                              '##jacquard.consensus.{0}DP_RANGE_ZSCORE.standard deviation_DP_range={1}\n'.format(JQ_CONSENSUS_TAG, rounded_std_dp)]
    meta_headers.extend(consensus_meta_headers)

    meta_headers.append(header)
    utils.write_output(writer, meta_headers, lines)
    logger.info("Wrote consensus-somatic-tagged VCF to [{}]", output_file)

# pylint: disable=C0301
def add_consensus(meta_headers, header, lines, writer, output_file):
    consensus_meta_headers = ['##FORMAT=<ID={0}SOM_SUM,Number=1,Type=Integer,Description="Jacquard consensus somatic call = sum(*{1}*)">\n'.format(JQ_CONSENSUS_TAG, utils.jq_somatic_tag), 
                              '##FORMAT=<ID={0}AF,Number=A,Type=Integer,Description="Jacquard consensus somatic call = average(*{1}*)">\n'.format(JQ_CONSENSUS_TAG, utils.jq_af_tag),
                              '##FORMAT=<ID={0}DP,Number=1,Type=Integer,Description="Jacquard consensus depth = average(*{1}*)">\n'.format(JQ_CONSENSUS_TAG, utils.jq_dp_tag)]

    meta_headers.extend(consensus_meta_headers)

    meta_headers.append(header)
    utils.write_output(writer, meta_headers, lines)

def roundTwoDigits(value): 
    new_values = []
    for val in value:
        if len(val.split(".")[1]) <= 2:
            new_values.append(val)
        else:
            new_values.append(str(round(100 * float(val))/100))
    return ",".join(new_values) 

def create_consensus_dict(key, val, input_dict, consensus_dict, type):
    split_som = val.split(",")
    for index, af in enumerate(split_som):
        if index in consensus_dict.keys():
            consensus_dict[index] += int(af) if type == "int" else float(af)
        else:
            consensus_dict[index] = int(af) if type == "int" else float(af)

    return consensus_dict

def get_consensus_som(field_dict):
    field_list = []
    for key in field_dict.keys():
        field_list.append(str(field_dict[key]))
    consensus = ",".join(field_list) if field_list != [] else 0

    return consensus

def get_consensus(consensus_tags, consensus_dict):
    if len(consensus_tags) != 0:
        average = []
        for key in consensus_dict.keys():
            avg = float(consensus_dict[key])/len(consensus_tags)
            rounded_avg = roundTwoDigits([str(avg)])
            average.append(str(rounded_avg))
        consensus = ",".join(average)
    else:
        consensus = 0.0

    return consensus

def get_range(consensus_tags, combined_dict, range):
    values = []
    if len(consensus_tags) > 1:
        for tag in consensus_tags:
            values.append(combined_dict[tag])
    if values != []:
        this_range = float(max(values)) - float(min(values))
        range.append(this_range)
    else:
        this_range = 0.0
        range.append(this_range)

    return range, this_range

def calculate_consensus(combined_dict, af_range, dp_range):
    consensus_af_tags = []
    consensus_dp_tags = []
    af = {}
    somatic = {}
    depth = {}

    for key in combined_dict.keys():
        if key.startswith("JQ_") and key.endswith(utils.jq_somatic_tag):
            if combined_dict[key] != ".":
                somatic = create_consensus_dict(key, combined_dict[key],
                                                combined_dict, somatic,
                                                "int")
        elif key.startswith("JQ_") and key.endswith(utils.jq_af_tag):
            if combined_dict[key] != ".":
                new_af = roundTwoDigits(combined_dict[key].split(","))
                af = create_consensus_dict(key, new_af, combined_dict, af,
                                           "float")
                consensus_af_tags.append(key)
        elif key.startswith("JQ_") and key.endswith(utils.jq_dp_tag):
            if combined_dict[key] != ".":
                depth = create_consensus_dict(key, combined_dict[key],
                                              combined_dict, depth,
                                              "float")
                consensus_dp_tags.append(key)

    consensus_af = get_consensus(consensus_af_tags, af)
    consensus_som = get_consensus_som(somatic)
    consensus_dp = get_consensus(consensus_dp_tags, depth)

    combined_dict[JQ_CONSENSUS_TAG + "AF"] = str(consensus_af)
    combined_dict[JQ_CONSENSUS_TAG + "SOM_SUM"] = str(consensus_som)
    combined_dict[JQ_CONSENSUS_TAG + "DP"] = str(consensus_dp)

    af_range, this_af_range = get_range(consensus_af_tags,
                                        combined_dict, af_range)
    dp_range, this_dp_range = get_range(consensus_dp_tags,
                                        combined_dict, dp_range)

    combined_dict[JQ_CONSENSUS_TAG + "AF_RANGE"] = str(this_af_range)
    combined_dict[JQ_CONSENSUS_TAG + "DP_RANGE"] = str(this_dp_range)

    return combined_dict, af_range, dp_range

def combine_format_values(format_col, sample, sample_column_name):
    new_format = [x + "_" + sample_column_name for x in format_col.split(":")]
    return OrderedDict(zip(new_format, sample.split(":")))

def process_line(line, af_range, af_mean, af_std, dp_range, dp_mean, dp_std, function):
    split_line = line.split("\t")
    format_col = split_line[8]
    samples = split_line[9:]
    new_samples = []

    for sample in samples:
        combined_dict = utils.combine_format_values(format_col, sample)

        if function == "zscore":
            combined_dict = calculate_zscore(af_mean, af_std, dp_mean,
                                             dp_std, combined_dict)
        elif function == "consensus":
            combined_dict, af_range, dp_range=calculate_consensus(combined_dict,
                                                                    af_range,
                                                                    dp_range)
        new_samples.append(":".join(combined_dict.values()))
    new_format = [":".join(combined_dict.keys())]
    new_line = "\t".join(split_line[:8] + new_format + new_samples) + "\n"

    return new_line

def add_subparser(subparser):
    parser = subparser.add_parser("consensus", help="Accepts a Jacquard-merged VCf file and creates a new file, adding consensus fields.")
    parser.add_argument("input", help="Path to Jacquard-merged VCF (or any VCF with Jacquard tags (e.g. JQ_SOM_MT)")
    parser.add_argument("output", help="Path to output VCf")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def execute(args, execution_context): 
    input_file = os.path.abspath(args.input)
    output_file = os.path.abspath(args.output)

    extension = os.path.splitext(os.path.basename(input_file))[1]
    if not os.path.isfile(input_file) or extension != ".vcf":
        logger.error("Input file [{}] must be a VCF file.", input_file)
        exit(1)

    extension = os.path.splitext(os.path.basename(output_file))[1]
    if extension != ".vcf":
        logger.error("Output file [{}] must be a VCF file.", output_file)
        exit(1)

    af_range = []
    dp_range = []

    tmp_file = output_file + ".tmp"
    input_file_reader = open(input_file, "r")
    tmp_file_writer = open(tmp_file, "w")

    logger.info("Adding consensus values to temporary file [{}]", tmp_file)
    meta_headers, header, lines = iterate_file(input_file_reader,
                                               af_range, dp_range, "consensus")
    add_consensus(meta_headers, header, lines, tmp_file_writer, tmp_file)

    input_file_reader.close()
    tmp_file_writer.close()

    tmp_file_reader = open(tmp_file, "r")
    output_file_writer = open(output_file, "w")

    logger.info("Adding z-scores to [{}]", output_file)

    meta_headers, header, lines = iterate_file(tmp_file_reader,
                                               af_range, dp_range, "zscore")
    add_zscore(meta_headers, header, lines, output_file_writer, output_file,
               af_range, dp_range)

    tmp_file_reader.close()
    output_file_writer.close()

    os.remove(tmp_file)
    logger.info("Removed temporary file [{}]", tmp_file)
