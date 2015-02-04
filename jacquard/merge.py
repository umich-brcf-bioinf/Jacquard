#pylint: disable=maybe-no-member, no-member, cell-var-from-loop, redefined-builtin
#pylint: disable=unnecessary-lambda
from __future__ import absolute_import
from collections import defaultdict, OrderedDict
import glob
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import *
import re

import jacquard.utils as utils
import jacquard.logger as logger

class PivotError(Exception):
    """Base class for exceptions in this module."""
    pass

class VariantPivoter():
    MISSING_REQUIRED_COLUMNS_ERROR="The columns of specified dataframe do " +\
        "not contain required columns {0}; review input data."

    def __init__(self, rows, combined_df=pd.DataFrame()):
        self._rows = rows
        self._combined_df = combined_df
        self.row_dict = {}

    def add_file(self, sample_file, header_index, caller, file_name=""):
        if file_name == "":
            file_name = sample_file

        initial_df  = create_initial_df(sample_file, header_index)

        if "#CHROM" in initial_df.columns:
            initial_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

        self.is_compatible(initial_df)
        fname_df, sample_columns = append_fname_to_samples(initial_df,
                                                           file_name,
                                                           self._rows,
                                                           caller)
        validate_sample_caller_vcfs(fname_df)

        self._combined_df = merge_samples(fname_df, self._combined_df,
                                          self._rows)

        return sample_columns

    def validate_sample_data(self):
        grouped = self._combined_df.groupby(self._rows)
        group = grouped.groups
#         for key, val in group.items():
#             if len(val) != 1:
#                 for column in self._combined_df:
#                     self.find_non_unique_rows(grouped, column, key, val)
        for column in self._combined_df:
            self.find_non_unique_cells(column)

        self._combined_df.reset_index(inplace=True)
        del self._combined_df["index"]

        return self._combined_df

    def find_non_unique_rows(self, grouped, column, key, val):
        col_data = []
        for index in val:
            data = grouped.get_group(key).ix[index, column]
            if data not in col_data:
                col_data.append(data)
        if len(col_data) != 1:
            self._combined_df.ix[val[-1], column] = "^"
            self._combined_df = self._combined_df.drop(self._combined_df.index[val[:-1]])

    def find_non_unique_cells(self, column):
        count = 0
        for data in self._combined_df[column]:
            if type(data) == np.ndarray:
                self._combined_df.ix[count, column] = "^"
            count += 1

    def sort_rows(self, df):
        df.reset_index(inplace=True)
        sorted_df = df.sort(self._rows)

        return sorted_df

    def is_compatible(self, initial_df):
        if "#CHROM" in initial_df.columns:
            initial_df.rename(columns={"#CHROM": "CHROM"}, inplace=True)

        initial_df = project_prepivot(initial_df)

        self._check_required_columns_present(initial_df)
        self._check_pivot_is_unique(initial_df)

        return initial_df

    def _check_required_columns_present(self, dataframe):
        required_columns = set(self._rows)

        if not required_columns.issubset(dataframe.columns.values):
            raise PivotError("Missing required columns; contact sysadmin.")

    def _check_pivot_is_unique(self, dataframe):
        group = self._rows
        grouped_df = dataframe.groupby(group)

        if len(grouped_df.groups) != len(dataframe):
            raise PivotError("Duplicate keys would result in an invalid pivot; "
                             "contact sysadmin.")

def validate_parameters(input_keys, first_line, header_names):
    invalid_fields = []
    fields = header_names.split("\t")

    for key in input_keys:
        if key not in fields:
            invalid_fields.append(key)

    message = "Invalid input parameter(s) "
    raise_err = 0
    if invalid_fields != []:
        message += str(invalid_fields)
        raise_err = 1

    return raise_err, message

def project_prepivot(df):
    remove_columns = set(["ID", "QUAL", "FILTER"])

    for col in df.columns:
        if col in remove_columns:
            del df[col]

    return df

def build_pivoter(sample_file, input_keys, header_index):
    initial_df  = create_initial_df(sample_file, header_index)

    pivoter = VariantPivoter(input_keys)
    pivoter.is_compatible(initial_df)

    return pivoter

def create_initial_df(sample_file, header_index):
    initial_df = pd.read_csv(sample_file, sep="\t", header=header_index,
                             dtype='str', mangle_dupe_cols=False)
    initial_df["INFO"] = "."

    return initial_df

def append_fname_to_samples(df, file_name, rows, caller):
    indexed_df = df.set_index(rows + ["INFO"])

    basename = os.path.basename(file_name)
    fname_prefix = basename.split(".")[0]
    new_df = indexed_df.rename(columns=lambda x: "|".join([caller,
                                                           fname_prefix, x]))
    sample_columns = new_df.columns.values
    reset_df = new_df.reset_index()

    return reset_df, sample_columns

def validate_sample_caller_vcfs(fname_df):
    columns = {}
    for col in fname_df.columns.values:
        if col in columns.keys():
            columns[col] += 1
        else:
            columns[col] = 1
    error = 0
    for key, val in columns.items():
        caller = key.split("|")[0]
        sample = "|".join(key.split("|")[1:])
        if val > 1:
            logger.error("Sample [{}] appears to be called by [{}] in "
                         "multiple files.", sample, caller)
            error = 1
    if error:
        raise utils.JQException("Some samples have calls for the same caller "
                                "in more than one file. Adjust or move problem "
                                "input files and try again.")

    return fname_df

def find_mult_alts(group):
    if len(group["ALT"]) > 1:
        group["INFO"] = "Mult_Alt"
    return group

def merge_samples(incoming_df, combined_df, rows):
    if combined_df.empty:
        combined_df = incoming_df
    else:
        combined_df = merge(combined_df, incoming_df, how="outer",
                            on=rows+["INFO"])

    return combined_df

def rearrange_columns(output_df):
    index = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    format_tags = []
    samples = []
    for i in list(output_df.columns.values):
        if i in index:
            continue
        elif re.search("FORMAT", i):
            format_tags.append(i)
        else:
            samples.append(i)

    headers = index + format_tags + samples

    ##change order of columns:
    output_df = output_df.ix[:,headers]
    return output_df

def create_dict(df, row, columns):
    #pylint: disable=anomalous-backslash-in-string
    file_dict = defaultdict(list)
    all_tags = []

    for column in columns:
        if re.search("\|", column) and not re.search("\|FORMAT", column):
            caller = column.split("|")[0]
            fname = column.split("|")[1]
            format_column = str(df.ix[row, "{0}|{1}|FORMAT".format(caller,
                                                                   fname)])
            sample_column = str(df.ix[row, column])

            format_column += ":sample_name"
            sample_column += ":" + column

#             format_sample = "{0}={1}".format(format_column, sample_column)
            format_sample_dict = _combine_format_values(format_column,
                                                             sample_column)
            tags = format_column.split(":")

#             key = "{0}|{1}".format(caller, fname)
#             file_dict[key].append(format_sample_dict)
            file_dict[caller].append(format_sample_dict)

            for tag in tags:
                try:
                    float(tag)
                except ValueError:
                    if tag not in all_tags:
                        all_tags.append(tag)

    return file_dict, all_tags

def add_all_tags(file_dict, all_sample_keys):
    for sample_list in file_dict.values():
        for sample_dict in sample_list:
            for samp_key in all_sample_keys:
                if samp_key not in sample_dict.keys():
                    sample_dict[samp_key] = "^"
                    #Keeps track of jacquard-introduced null values.

    return file_dict

def sort_format_tags(file_dict):
    sorted_file_dict = defaultdict(list)
    for fname, sample_dicts in file_dict.items():
        new_samp_dicts = []
        for sample_dict in sample_dicts:
            new_sample_dict = OrderedDict(sorted(sample_dict.iteritems()))
            new_samp_dicts.append(new_sample_dict)
        sorted_file_dict[fname] = new_samp_dicts

    return sorted_file_dict

def remove_non_jq_tags(df, file_dict):
    #pylint: disable=unused-argument
    sample_keys = []
#     sample_names = {}
    for sample_list in file_dict.values():
        for sample in sample_list:
            for key in sample.keys():
#                 sample_names[sample["sample_name"]] = 1
                if not re.search("JQ_", key) and key != "sample_name":
                    del sample[key]
                else:
                    if key not in sample_keys:
                        sample_keys.append(key)
    file_dict = add_all_tags(file_dict, sample_keys)
    sort_format_tags(file_dict)

    return file_dict

def cleanup_df(df, file_dict):
    for key in file_dict.keys():
        try:
            del df[key + "|FORMAT"]
        except KeyError:
            pass
    for col in df.columns:
        df = df.applymap(lambda x: str(x).replace(":" + col, ""))
        df = df.applymap(lambda x: str(x).replace(col + ":", ""))
        df = df.applymap(lambda x: str(x).replace(col, "."))

    df = df.applymap(lambda x: str(x).replace(":sample_name", ""))
    df = df.applymap(lambda x: str(x).replace("sample_name:", ""))
    df = df.applymap(lambda x: str(x).replace("sample_name", "."))

    df.replace("nan", ".", inplace=True)

    return df

def create_merging_dict(df, row, columns):
    #pylint: disable=anomalous-backslash-in-string
    file_dict = defaultdict(list)
    sample_names = []
    sample_columns = []

    format_column = str(df.ix[row, "FORMAT"])

    for column in columns:
        if re.search("\|", column) and not re.search("\|FORMAT", column):
            fname = column.split("|")[1]
            sample = column.split("|")[2]
            sample_column = str(df.ix[row, column])
            sample_columns.append(sample_column)
            sample_names.append(column)

            format_sample_dict = _combine_format_values(format_column,
                                                             sample_column)

            file_key = "{0}|{1}".format(fname, sample)
            if format_sample_dict != OrderedDict():
                file_dict[file_key].append(format_sample_dict)

    return file_dict

def remove_old_columns(df):
    #pylint: disable=anomalous-backslash-in-string
    columns_to_remove = []
    for row, col in df.T.iteritems():
        columns = col.index.values
        for column in columns:
            if re.search("\|", column) and len(column.split("|")) == 3:
                if column not in columns_to_remove:
                    columns_to_remove.append(column)

    for col in columns_to_remove:
        del df[col]

    return df

def combine_format_columns(df, all_inconsistent_sample_sets):
    #pylint: disable=unused-argument
    row_total = len(df.index)
    logger.info("Processing merged matrix phase 1: [{} x {}] rows x columns",
                row_total, len(df.columns))
    for row, col in df.T.iteritems():
        columns = col.index.values
        if row % 1000 == 1:
            logger.info("Processing merged matrix phase 1: [{}/{}] rows "
                        "processed ({}% complete)", row, row_total,
                        int(100 * row/row_total))
        file_dict, all_tags = create_dict(df, row, columns)
        file_dict = remove_non_jq_tags(df, file_dict)

        for key, val in file_dict.items():
            for samp_dict in val:
                format_tags = []
                sample = []
                sorted_dict = OrderedDict(sorted(samp_dict.items()))

                for sorted_key, sorted_val in sorted_dict.items():
                    format_tags.append(sorted_key)
                    sample.append(sorted_val)

                df.ix[row, "FORMAT"] = ":".join(format_tags)
                df.ix[row, samp_dict["sample_name"]] = ":".join(sample)

    df = cleanup_df(df, file_dict)
    logger.info("Processing merged matrix: phase 2 [{}] rows", len(df.index))
    for row, col in df.T.iteritems():
        if row % 1000 == 1:
            logger.info("Processing merged matrix phase 2: [{}/{}] rows "
                        "processed ({}%)", row, row_total,
                        int(100 * row/row_total))
        columns = col.index.values
        file_dict = create_merging_dict(df, row, columns)
#         merge_by = len(file_dict.items())
        for key, vals in file_dict.items():
            complete_dict = {}

            for val_dict in vals:
                for format_key, format_val in val_dict.items():
                    if format_val == "^":
                        format_val = "."
                    if format_key in complete_dict:
                        if complete_dict[format_key]==".":
                            complete_dict[format_key] = format_val
                    else:
                        complete_dict[format_key] = format_val
            sorted_complete_dict = OrderedDict(sorted(complete_dict.items()))

            df.ix[row, key] = ":".join(sorted_complete_dict.values())

    df = remove_old_columns(df)

    return df

def determine_merge_execution_context(all_merge_context, all_merge_column_context, sample_columns, sample_file, count):
    samples = []
    actual_sample_columns = []

    samp_count = 0
    for samp_column in sample_columns:
        samp_name = samp_column.split("|")[-1]
        if samp_name != "FORMAT":
            samp_count += 1
            samples.append(samp_name)
            actual_sample_columns.append(samp_column)
            merge_column_context = "##jacquard.merge.sample_column{0}={1}({2})".format(samp_count, samp_column, os.path.basename(sample_file))
            all_merge_column_context.append(merge_column_context)

    merge_context = "##jacquard.merge.file{0}={1}({2})".format(count, os.path.basename(sample_file), samples)
    all_merge_context.append(merge_context)

    return all_merge_context, all_merge_column_context

def print_new_execution_context(execution_context, out_file):
    execution_context[-1] = execution_context[-1] + "\n"

    out_file.write("\n".join(execution_context))
    out_file.close()

def create_new_line(alt_allele_number, fields):
    alt = fields[4].split(",")[alt_allele_number]
    format_tags = fields[8]
    samples = fields[9:]
    new_samples = []

    for sample in samples:
        format_sample_dict = _combine_format_values(format_tags, sample)
        new_dict = OrderedDict()
        for key, val in format_sample_dict.items():
            #only care about splitting jacquard tags
            if re.search("JQ_.*AF", key):
                split_val = val.split(",")
                if len(split_val) > 1:
                    new_dict[key] = split_val[alt_allele_number]
                else:
                    new_dict[key] = val
            else:
                new_dict[key] = val
        new_samples.append(":".join(new_dict.values()))

    new_line = fields[0:4] + [alt] + fields[5:9] + new_samples

    return "\t".join(new_line) + "\n"

def determine_caller_and_split_mult_alts(reader, writer, unknown_callers):
    caller = "unknown"
    for line in reader:
        if line.startswith("##jacquard.tag.caller="):
            caller = line.split("=")[1].rstrip()
            writer.write(line)
        elif line.startswith("#"):
            writer.write(line)
        else:
            fields = line.split("\t")
            alt = fields[4]
            alts = alt.split(",")
            if len(alts) > 1: #there's a mult-alt
                count = 0
                for _ in alts:
                    new_line = create_new_line(count, fields)
                    writer.write(new_line)
                    count += 1
            else:
#                 new_line = fields[:7] + ["."] + fields[8:]
#                 writer.write("\t".join(new_line))
                writer.write("\t".join(fields))

    if caller == "unknown":
        logger.error("Unable to determine variant caller for file [{}]", reader)
        unknown_callers += 1

    return caller, unknown_callers

def validate_samples_for_callers(all_merge_column_context, all_inconsistent_sample_sets):
    sample_dict = defaultdict(list)
    samples = []
    for message in all_merge_column_context:
        message_info = message.split("=")[1]
        column = message_info.split("(")[0]
        fname = message_info.split("(")[1].strip(")")
        caller = column.split("|")[0]
        sample = column.split("|")[1]
        sample_column = column.split("|")[2]

        samples.append(sample)
        sample_dict[caller].append(sample)
    logger.info("Detected VCFs from {}", sample_dict.keys())

    warn = 0
    for key, val in sample_dict.items():
        missing = []
        for sample in samples:
            if sample not in val:
                missing.append(sample)
        if missing != []:
            logger.warning("Samples {} were not called by {}", missing, key)
            warn = 1

    if warn == 1 and all_inconsistent_sample_sets == False:
        raise utils.JQException("Some samples were not present for all "
                                "callers. Review log warnings and move/adjust "
                                "input files as appropriate.")
    elif warn == 1 and all_inconsistent_sample_sets == True:
        logger.warning("Some samples were not present for all callers.")

    return 1

def _add_mult_alt_flags(df):
    return df.groupby(by=["CHROM", "POS", "REF"]).apply(find_mult_alts)

def process_files(sample_file_readers, input_file, output_path, input_keys, headers, header_names, first_line, all_inconsistent_sample_sets, execution_context, pivot_builder=build_pivoter):
    first_file_reader = sample_file_readers[0]
    first_file      = first_file_reader

    raise_err, message = validate_parameters(input_keys,
                                             first_line,
                                             header_names)
    if raise_err == 1:
        raise PivotError(message)

    pivoter  = pivot_builder(first_file, input_keys, headers[0])

    all_merge_context = []
    all_merge_column_context = []
    unknown_callers = 0

    output_file = os.path.dirname(output_path)
    new_dir = os.path.join(output_file, "splitMultAlts")
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    logger.info("Splitting mult-alts in input files. Using [{}] as input "
                "directory.", new_dir)

    total_number_of_files = len(sample_file_readers)
    count = 1
    for sample_file in sample_file_readers:
        logger.info("Reading [{}] ({}/{})", os.path.basename(sample_file),
                    count, total_number_of_files)
        fname, extension = os.path.splitext(os.path.basename(sample_file))

        new_sample_file = os.path.join(new_dir,
                                       fname + ".splitMultAlts" + extension)

        reader = open(sample_file, "r")
        writer = open(new_sample_file, "w")
        caller, unknown_callers = determine_caller_and_split_mult_alts(reader,
                                                                       writer,
                                                                       unknown_callers)
        reader.close()
        writer.close()

        sample_columns = pivoter.add_file(new_sample_file, headers[count - 1],
                                          caller)
        count += 1

        all_merge_context, all_merge_column_context = determine_merge_execution_context(all_merge_context, all_merge_column_context, sample_columns, new_sample_file, count)
    if unknown_callers != 0:
        raise utils.JQException("Unable to determine variant caller for [{}] input files. Run (jacquard tag) first.", unknown_callers)

#     pivoter._combined_df = merge(pivoter._combined_df, pd.DataFrame(columns=pivoter._rows+["INFO"]), how="outer", on=pivoter._rows+["INFO"])
#     print pivoter._combined_df

    pivoter._combined_df = _add_mult_alt_flags(pivoter._combined_df)

    validate_samples_for_callers(all_merge_column_context,
                                 all_inconsistent_sample_sets)

    new_execution_context = all_merge_context + all_merge_column_context

    execution_context.extend(new_execution_context)
#     writer = open(output_path, "w")
    logger.info("Merging sample data (this may take a while)")

    logger.info("Merging sample data: validating sample data (1/6)")
    pivoter.validate_sample_data()

    logger.info("Merging sample data: combining format columns (2/6)")
    formatted_df = combine_format_columns(pivoter._combined_df, all_inconsistent_sample_sets)

    logger.info("Merging sample data: rearranging columns (3/6)")
    rearranged_df = rearrange_columns(formatted_df)
    try:
        rearranged_df.ix[:, "CHROM"] = rearranged_df.ix[:, "CHROM"].apply(lambda x: int(x.strip("chr")))
    except ValueError:
        rearranged_df.ix[:, "CHROM"] = rearranged_df.ix[:, "CHROM"].apply(lambda x: x.strip("chr"))
    rearranged_df.ix[:, "POS"] = rearranged_df.ix[:, "POS"].apply(lambda x: int(x))

    logger.info("Merging sample data: sorting rows (4/6)")
    sorted_df = pivoter.sort_rows(rearranged_df)
    sorted_df.ix[:, "CHROM"] = sorted_df.ix[:, "CHROM"].apply(lambda x: "chr" + str(x))
    sorted_df.ix[:, "POS"] = sorted_df.ix[:, "POS"].apply(lambda x: str(x))
    
    if "index" in sorted_df:
        del sorted_df["index"]

    logger.info("Merging sample data: cleanup (5/6)")
    sorted_df = sorted_df.fillna(".")
    sorted_df.rename(columns={"CHROM": "#CHROM"}, inplace=True)

    logger.info("Merging sample data: saving (6/6)")
    with open(output_path, "a") as f:
        f.write("\n".join(execution_context)+"\n")
        sorted_df.to_csv(f, index=False, sep="\t")

    logger.info("Merged [{}] VCf files to [{}]",
                len(sample_file_readers),
                output_path)

def determine_input_keys(input_file):
    for file_name in listdir(input_file):
        if isfile(join(input_file, file_name)):
            fname, extension = os.path.splitext(file_name)
            if extension == ".vcf":
                return ["CHROM", "POS", "REF", "ALT"]

            else:
                raise PivotError("Cannot determine columns to be used as keys "
                                 "for the pivoting from file [{0}]. Please "
                                 "specify parameter [-k] or "
                                 "[--keys]".format(os.path.abspath(file_name)))
    
def get_headers_and_readers(in_files):
    sample_file_readers = []
    headers = []
    header_names = []
    first_line = []
    meta_headers = []
    invalid_files = []
    for in_file in in_files:
        f = open(in_file, 'r')
        count = -1
        invalid = 1
        for line in f:
            count += 1
            if line.startswith("##"):
                if re.search("##FORMAT=<ID=JQ_", line):
                    meta_headers.append(line.rstrip())
                if re.search("##jacquard.tag.caller=", line):
                    invalid = 0
            elif line.startswith("#"):
                headers.append(count)
                header_names.append(line)
            else:
                first_line.append(line)
                break
        if invalid:
            invalid_files.append(in_file)

        f.close()
        sample_file_readers.append(in_file)

    if invalid_files:
        raise utils.JQException("VCF file(s) [{}] have no Jacquard tags. Run "
                                "[jacquard tag] on these files and try again.",
                                invalid_files)

    header_names = header_names[0]
    header_names = re.sub(r'#CHROM', 'CHROM', header_names)

    return sample_file_readers, headers, header_names, first_line, meta_headers

def add_subparser(subparser):
    parser = subparser.add_parser("merge", help="Accepts a directory of VCFs and returns a single merged VCF file.")
    parser.add_argument("input", help="Path to directory containing VCFs. Other file types ignored")
    parser.add_argument("output", help="Path to output variant-level VCF file")
    parser.add_argument("-k", "--keys",
        help="Columns to be used as keys for the pivoting. Default keys for VCF are CHROM,POS,ID,REF,ALT,QUAL,FILTER.")
    parser.add_argument("-a", "--allow_inconsistent_sample_sets", action="store_true", default=False, help="Allow inconsistent sample sets across callers. Not recommended.")
    parser.add_argument("-v", "--verbose", action='store_true')
    parser.add_argument("--force", action='store_true', help="Overwrite contents of output directory")

def get_required_input_output_types():
    return ("foo","bar")

def report_prediction(args):
    return "foo"

def execute(args, execution_context):
    input_file = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)

    output_file, outfile_name = os.path.split(output_path)
    utils.validate_directories(input_file, output_file)
    input_keys = args.keys.split(",") if args.keys else determine_input_keys(input_file)
    all_inconsistent_sample_sets = args.allow_inconsistent_sample_sets

    fname, extension = os.path.splitext(outfile_name)
#     if extension != ".vcf": 
#         raise utils.JQException("Error. Specified output {} must have a .vcf extension", output_path)

    in_files = sorted(glob.glob(os.path.join(input_file,"*.vcf")))
    if len(in_files) < 1:
        raise utils.JQException("Error: Specified input directory [{}] "
                                "contains no VCF files. Check parameters and "
                                "try again.")

    sample_file_readers, headers, header_names, first_line, meta_headers = get_headers_and_readers(in_files)

    execution_context.extend(meta_headers + ["##fileformat=VCFv4.2"])

    process_files(sample_file_readers, input_file, output_path, input_keys,
                  headers, header_names, first_line,
                  all_inconsistent_sample_sets, execution_context)

def _combine_format_values(format_field, sample):
    return OrderedDict(zip(format_field.split(":"), sample.strip().split(":")))

