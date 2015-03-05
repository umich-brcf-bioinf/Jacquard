from openpyxl import load_workbook
from openpyxl.style import Color, Fill, Font, Alignment
import os
import re

class PivotError(Exception):
    """Base class for exceptions in this module."""
    pass

def style_workbook(input_file, output_file):
    workbook = load_workbook(input_file)
    worksheet = workbook.active

    for col in worksheet.columns: #pylint: disable=no-member
        if col[0].value == "IGV":
            format_links(worksheet, col, "IGV")
        elif col[0].value == "UCSC":
            format_links(worksheet, col, "UCSC")
        elif col[0].value == "dbSNP":
            format_links(worksheet, col, "dbSNP")
#         count = 0
#         colors = ["FFFF99", "FFCC00", "FF9900", "FF6600"]
        # for tag in pivot_values:
            # if re.search(tag, col[0].value):
                # fill_cell(col, colors[count])
            # count += 1

        annotation_headers = ("CHROM|POS|ID|REF|ALT|Mult_Alt|Mult_Gene|"
                              "ANNOTATED_ALLELE|GENE_SYMBOL|IGV|UCSC|dbSNP|"
                              "QUAL|FILTER|INFO")
        if re.search(annotation_headers, col[0].value):
            fill_cell(col, "C0C0C0")

        if col[0].style.fill.start_color.index == "FFFFFFFF":
            fill_cell(col, "D7E8F0")
        if col[0].value == "Mult_Alt":
            fill_mults(worksheet, col)
        if col[0].value == "Mult_Gene":
            fill_mults(worksheet, col)
        if col[0].value == "INFO":
            fill_mult_alts(worksheet, col)

    workbook.save(output_file)

def format_links(worksheet, column, cell_value):
    for pos in column:
        if pos != column[0]:
            desired_cell = str(pos).strip("(<>)").split(".")[-1]
            worksheet[desired_cell].hyperlink = pos.value
            pos.value = cell_value
            pos.style.font.color.index = Color.BLUE
            pos.style.font.underline = Font. UNDERLINE_SINGLE

def fill_cell(column, color):
    column[0].style.alignment.text_rotation = 90
    column[0].style.alignment.vertical = Alignment.VERTICAL_BOTTOM
    column[0].style.fill.fill_type = Fill.FILL_SOLID
    column[0].style.fill.start_color.index = color

def fill_row(row, color):
    row.style.fill.fill_type = Fill.FILL_SOLID
    row.style.fill.start_color.index = color

def fill_mults(worksheet, col):
    for row in col:
        if row.value == "True":
            coordinate = row.address
            row = row.address[1:]
            for item in worksheet.range("A" + row + ":" + coordinate):
                for cell in item:
                    if cell.value != "None":
                        fill_row(cell, "FFD698")

def fill_mult_alts(worksheet, col):
    for row in col:
        if row.value == "Mult_Alt":
            #coordinate = row.address
            row = row.address[1:]
#             for item in worksheet.range("A" + row + ":" + coordinate):
            for item in worksheet.range("A" + row + ":D" + row):
                for cell in item:
                    if cell.value != "None":
                        fill_row(cell, "FFD698")

def add_subparser(subparser):
    #pylint: disable=line-too-long
    parser_pivot = subparser.add_parser("style", help="Accepts an XLSX file and outputs a styled XLSX file with conditional formatting")
    parser_pivot.add_argument("input_file", help="Input XLSX file that has been formatted by Jacquard")
    parser_pivot.add_argument("output_file", help="Output XLSX file")

def execute(args, dummy):
    input_file = os.path.abspath(args.input_file)
    output_file = os.path.abspath(args.output_file)

    for file_name in (input_file, output_file):
        dummy, extension = os.path.splitext(file_name)
        if extension != ".xlsx":
            print ("Error. Specified output {0} must "
                   "have a .vcf extension").format(file_name)
            exit(1)

    if not os.path.isfile(input_file):
        print ("Error. Specified input directory {0} "
               "does not exist").format(input_file)
        exit(1)

#    script_dir = os.path.dirname(os.path.abspath(__file__))
#     pd.set_option('chained_assignment', None)

    style_workbook(input_file, output_file)
