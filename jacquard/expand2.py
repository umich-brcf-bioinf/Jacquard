import os

import vcf as vcf
import utils as utils

def _get_vcf_reader(input_file):
    file_reader = vcf.FileReader(input_file)
    vcf_reader = vcf.VcfReader(file_reader)
    
    return vcf_reader

def _parse_meta_headers(meta_headers):
    info_fields = []
    format_tags = []
    for meta_header in meta_headers:
        if meta_header.startswith("##FORMAT="):
            tag = meta_header.split(",")[0].split("=")[-1]
            format_tags.append(tag)
        elif meta_header.startswith("##INFO="):
            field = meta_header.split(",")[0].split("=")[-1]
            info_fields.append(field)
            
    info_fields.sort()
    format_tags.sort()
    
    if len(info_fields) == 0 or len(format_tags) == 0:
        raise utils.JQException("Unable to parse meta_headers for INFO and/or " + 
                                "FORMAT fields. Review input and try again.")
    
    return info_fields, format_tags

def _append_format_tags_to_samples(format_tags, samples):
    format_sample_header = []
    
    for sample in samples:
        for tag in format_tags:
            modified_field = tag + "|" + sample
            format_sample_header.append(modified_field)
            
    return format_sample_header

def _get_headers(vcf_reader):
    split_column_header = vcf_reader.column_header.split("\t")
    
    column_header_no_samples = split_column_header[0:7]
    samples = split_column_header[9:]
    
    (info_header, format_tags) = _parse_meta_headers(vcf_reader.metaheaders)
    format_sample_header = _append_format_tags_to_samples(format_tags, samples)
    
    return column_header_no_samples, info_header, format_sample_header

def _parse_info_field(vcf_record, info_header):
    info_dict = vcf_record.get_info_dict()
    write_info_columns = []
    for tag in info_header:
        info_cell = ""
        if tag in info_dict:
            info_cell = info_dict[tag]
        write_info_columns.append(info_cell)
    return write_info_columns
    
def _write_vcf_records(vcf_reader, file_writer, info_header):
    for record in vcf_reader.vcf_records():
        write_info_columns = _parse_info_field(record, info_header)
        write_row = [record.chrom,record.pos,record.id,record.ref,record.alt,record.qual,record.filter]
        write_row.extend(write_info_columns)
        file_writer.write("\t".join(write_row))

def add_subparser(subparser):
    parser_pivot = subparser.add_parser("expand2", help="Pivots annotated VCF file so that given sample specific information is fielded out into separate columns. Returns an Excel file containing concatenation of all input files.")
    parser_pivot.add_argument("input_file", help="Path to annotated VCF. Other file types ignored")
    parser_pivot.add_argument("output_file", help="Path to output variant-level XLSX file")
    parser_pivot.add_argument("-v", "--verbose", action='store_true')
        
def execute(args, execution_context):
    input_file = os.path.abspath(args.input_file)
    output_path = os.path.abspath(args.output_file)
    
    output_dir = os.path.split(output_path)[0]
    try:
        os.mkdir(output_dir)
    except:
        pass  
    
    vcf_reader = _get_vcf_reader(input_file)
    column_header_no_samples, info_header, format_sample_header = _get_headers(vcf_reader)
    
    complete_header = column_header_no_samples + info_header + format_sample_header

    file_writer = vcf.FileWriter(output_path)
    file_writer.write(complete_header)
    
    _write_vcf_records(vcf_reader, file_writer, info_header)
    
    file_writer.close()
    
    