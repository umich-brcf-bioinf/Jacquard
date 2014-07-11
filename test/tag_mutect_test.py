#!/usr/bin/python2.7
import os
import unittest
from bin.tag_mutect import AlleleFreqTag, DepthTag, SomaticTag, LineProcessor, FileProcessor

class AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_AF_MT,Number=1,Type=Float, Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA).">', AlleleFreqTag().metaheader)
                
    def test_format_missingAFTag(self):
        tag = AlleleFreqTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format_rounds(self):
        tag = AlleleFreqTag()
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.2:0.2"), tag.format("A:FA", "1:0.2"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.20:0.20"), tag.format("A:FA", "1:0.20"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.204:0.2"), tag.format("A:FA", "1:0.204"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.205:0.21"), tag.format("A:FA", "1:0.205"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:0.206:0.21"), tag.format("A:FA", "1:0.206"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:1.0:1.0"), tag.format("A:FA", "1:1.0"))
        self.assertEqual(("A:FA:JQ_AF_MT", "1:1.00:1.00"), tag.format("A:FA", "1:1.00"))

class DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_DP_MT,Number=1,Type=Float, Description="Jacquard depth for MuTect (based on DP).">', DepthTag().metaheader)
                
    def test_format_missingDPTag(self):
        tag = DepthTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format(self):
        tag = DepthTag()
        self.assertEqual(("A:DP:JQ_DP_MT", "1:42:42"), tag.format("A:DP", "1:42"))

class SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID=JQ_SOM_MT,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1= somatic (based on SS FORMAT tag).">', SomaticTag().metaheader)
                
    def test_format_missingSSTag(self):
        tag = SomaticTag()
        format_param_string = "A:B"
        format_value_string = "1:2"
        self.assertEqual(("A:B", "1:2"), tag.format(format_param_string, format_value_string))
                
    def test_format(self):
        tag = SomaticTag()
        self.assertEqual(("A:SS:JQ_SOM_MT", "1:2:1"), tag.format("A:SS", "1:2"))
        self.assertEqual(("A:SS:JQ_SOM_MT", "1:1:0"), tag.format("A:SS", "1:1"))

class LineProcessorTestCase(unittest.TestCase):
    def test_process_line_singleSample(self):
        tag = MockLowerTag()    
        processor = LineProcessor([tag])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|X:Y:Z".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        actual_format_params = actual_line.split("\t")[8]
        actual_format_values = actual_line.split("\t")[9:]
        
        self.assertEqual("a:b:c", actual_format_params)
        self.assertEqual(["x:y:z"], actual_format_values)
        self.assertEqual("chr1\t42\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\tx:y:z", actual_line)
        
    def test_process_line_MultSample(self):
        tag = MockLowerTag()    
        processor = LineProcessor([tag])
        input_line = "chr1|42|.|ref|alt|qual|filter|INFO|A:B:C|U:V:W|X:Y:Z".replace("|", "\t")

        actual_line = processor.add_tags(input_line)
        actual_format_params = actual_line.split("\t")[8]
        actual_format_values = actual_line.split("\t")[9:]
        
        self.assertEqual("a:b:c", actual_format_params)
        self.assertEqual(["u:v:w", "x:y:z"], actual_format_values)
        self.assertEqual("chr1\t42\t.\tref\talt\tqual\tfilter\tINFO\ta:b:c\tu:v:w\tx:y:z", actual_line)

class FileProcessorTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('''##jacquard.version=0.1\n
        ##jacquard.tagMutect.command=tagMutect inDir outDir\n
        ##jacquard.tagMutect.cwd=/foo/bar''', FileProcessor().metaheader)
         
#     def test_process_file(self):
#         script_dir = os.path.dirname(os.path.abspath(__file__))
#         input_file = script_dir + "/test_input/P2_tiny_test_input.txt"
#         processor = FileProcessor(input_file)
# 
#         actual_line = processor.get_line()
#         expected_line = "1       1       17730   .       C       A       A       .       .       WASH7P  .       .       .       .       HOM     .       SNP     Novel_locus     .       .         .       HIGH    SPLICE_SITE_ACCEPTOR    .       SPLICE_SITE_ACCEPTOR(HIGH|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000430492|7|1),EXON(MODIFIER|||||WAS         H7P|unprocessed_pseudogene|NON_CODING|ENST00000423562|4|1),EXON(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000438504|5|1),EXON(MODIFIER|||||WASH         7P|unprocessed_pseudogene|NON_CODING|ENST00000488147|5|1),EXON(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000537342|4|1),EXON(MODIFIER|||||WASH7         P|unprocessed_pseudogene|NON_CODING|ENST00000538476|5|1),EXON(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000541675|4|1),INTRON(MODIFIER|||||WASH         7P|unprocessed_pseudogene|NON_CODING|ENST00000430492|6|1)       POPULATION_SNPFreq=Novel_locus;POPULATION_AF_CATEGORY=Novel_locus;POPULATION_AF=.;POPULATION_AF_S         OURCE=.;NUCMATCH=0;in_strand=2;in_rna_coverage=11;in_rna_meanq=37.36;in_rna_base_count=2,9,0,0;in_rna_all_subs=CA;in_rna_frequency=0.18;in_dna_coverage=292;in_dn         a_meanq=37.09;in_dna_base_count=8,284,0,0;in_dna_all_subs=-;in_dna_frequency=0.00;in_rdd=1;SNP;HOM;VARTYPE=SNP;dbNSFP_rollup_tolerated=0;dbNSFP_rollup_damaging=0         ;dbNSFP_rollup_total=0  .       .       .       .       .       .       .       .       .       .       .       0       0       0       RDD_GT:RNA_COVERAGE:RNA_F         REQ:RNA_ACGT    C_A:11:0.18:2,9,0,0     .       .\n"
#         
#         self.assertEqual(expected_line, actual_line)

        
class MockLowerTag():
    def format(self, params, values):
        return (params.lower(), values.lower())
    