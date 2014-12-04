import unittest

import jacquard.variant_callers.mutect as mutect
from jacquard.utils import __version__, JQException
from jacquard.vcf import VcfRecord 

class MockWriter():
    def __init__(self):
        self._content = []
        self.opened = False
        self.closed = False

    def open (self):
        self.opened = True

    def write(self, content):
        self._content.extend(content.splitlines())
        
    def lines(self):
        return self._content

    def close(self):
        self.closed = True

class MockReader():
    def __init__(self, lines = []):
        self._lines_iter = iter(lines)
        self.opened = False
        self.closed = False
        self.input_filepath = "foo"

    def open (self):
        self.opened = True

    def read_lines(self):
        return self._lines_iter

    def close(self):
        self.closed = True


class AlleleFreqTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}AF,Number=A,Type=Float,Description="Jacquard allele frequency for MuTect: Decimal allele frequency rounded to 2 digits (based on FA)",Source="Jacquard",Version={1}>'.format(mutect.JQ_MUTECT_TAG, __version__), mutect._AlleleFreqTag().metaheader)

    def test_format_missingAFTag(self):
        tag = mutect._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        originalVcfRecord = VcfRecord(line)
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())

    def test_format_presentAFTag(self):
        tag = mutect._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3|0.567:SA.2:SA.3|0.834:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3:{0}AF|0.567:SA.2:SA.3:0.57|0.834:SB.2:SB.3:0.83\n".format(mutect.JQ_MUTECT_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_multAlt(self):
        tag = mutect._AlleleFreqTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3|0.5,0.8:SA.2:SA.3|0.7,0.6:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FA:F2:F3:{0}AF|0.5,0.8:SA.2:SA.3:0.5,0.8|0.7,0.6:SB.2:SB.3:0.7,0.6\n".format(mutect.JQ_MUTECT_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class DepthTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}DP,Number=1,Type=Float,Description="Jacquard depth for MuTect (based on DP)",Source="Jacquard",Version={1}>'.format(mutect.JQ_MUTECT_TAG, __version__), mutect._DepthTag().metaheader)

    def test_format_missingDPTag(self):
        tag = mutect._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        originalVcfRecord = VcfRecord(line)
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(originalVcfRecord.asText(), processedVcfRecord.asText())

    def test_format_presentDPTag(self):
        tag = mutect._DepthTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3|2:SA.2:SA.3|4:SB.2:SB.3\n".replace('|',"\t")
        expected = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|DP:F2:F3:{0}DP|2:SA.2:SA.3:2|4:SB.2:SB.3:4\n".format(mutect.JQ_MUTECT_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class SomaticTagTestCase(unittest.TestCase):
    def test_metaheader(self):
        self.assertEqual('##FORMAT=<ID={0}HC_SOM,Number=1,Type=Integer,Description="Jacquard somatic status for MuTect: 0=non-somatic,1=somatic (based on SS FORMAT tag)",Source="Jacquard",Version={1}>'.format(mutect.JQ_MUTECT_TAG, __version__), mutect._SomaticTag().metaheader)

    def test_format_missingSSTag(self):
        tag = mutect._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:{0}HC_SOM|SA.1:SA.2:SA.3:0|SB.1:SB.2:SB.3:0\n").format(mutect.JQ_MUTECT_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

    def test_format_presentSSTag(self):
        tag = mutect._SomaticTag()
        line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3|2:SA.2:SA.3|5:SB.2:SB.3\n".replace('|',"\t")
        expected = ("CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|SS:F2:F3:{0}HC_SOM|2:SA.2:SA.3:1|5:SB.2:SB.3:0\n").format(mutect.JQ_MUTECT_TAG).replace('|',"\t")
        processedVcfRecord = VcfRecord(line)
        tag.format(processedVcfRecord)
        self.assertEquals(expected, processedVcfRecord.asText())

class MockTag(object):
    def __init__(self, field_name, field_value, metaheader=None):
        self.field_name = field_name
        self.field_value = field_value
        self.metaheader = metaheader

    def format(self, vcfRecord):
        vcfRecord.insert_format_field(self.field_name, {0:self.field_value, 1:self.field_value})

class Mutect_TestCase(unittest.TestCase):
    def setUp(self):
        self.caller = mutect.Mutect()

    def test_validate_vcfs_in_directory(self):
        in_files = ["A.vcf","B.vcf"]
        self.caller.validate_vcfs_in_directory(in_files)

        in_files = ["A.vcf","B"]
        self.assertRaisesRegexp(JQException, "ERROR: Non-VCF file in directory. Check parameters and try again", self.caller.validate_vcfs_in_directory, in_files)

    def test_decorate_files(self):
        filenames = ["A/A.vcf"]
        decorator = "normalized"
        actual_filenames = self.caller.decorate_files(filenames, decorator)
        expected_filenames = "A.normalized.vcf"
        self.assertEquals(expected_filenames,actual_filenames)

    def test_validateInputFile_isValid(self):
        metaheaders = ["##MuTect=blah"]
        self.assertTrue(self.caller.validate_input_file(metaheaders, "#column_header"))

    def test_validateInputFile_isNotValid(self):
        metaheaders = ["Foo"]
        self.assertFalse(self.caller.validate_input_file(metaheaders, "#column_header"))

    def test_addTags(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        self.caller.tags=[MockTag("mockTag", 42)]
        actual_line = self.caller.add_tags(VcfRecord(input_line))

        expected_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3:mockTag|SA.1:SA.2:SA.3:42|SB.1:SB.2:SB.3:42\n".replace('|',"\t")

        self.assertEquals(expected_line, actual_line)

    def test_updateMetaheader(self):
        self.caller.tags=[MockTag("mockTag", 42, "##my_metaheader\n")]
        actual_metaheader = self.caller.get_new_metaheaders()

        self.assertEquals(["##my_metaheader\n"], actual_metaheader)

    def test_normalize(self):
        writer = MockWriter()
        content = ["foo", "bar", "baz"]
        reader = MockReader(content)
        self.caller.normalize(writer,[reader])

        self.assertTrue(reader.opened)
        self.assertTrue(reader.closed)
        self.assertTrue(writer.opened)
        self.assertTrue(writer.closed)
        self.assertEquals(content, writer.lines())

    def test_normalize_raisesExceptionIfTwoInputFiles(self):
        self.assertRaisesRegexp(JQException, r"MuTect .* but found \[2\]\.", self.caller.normalize, MockWriter(), [MockReader(), MockReader()])

    def test_normalize_changes_column_headers(self):
        writer = MockWriter()
        content = ["##MuTect=foo normal_sample_name=normal_sample tumor_sample_name=tumor_sample foo=bar", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor_sample\tnormal_sample"]
        reader = MockReader(content)
        self.caller.normalize(writer,[reader])

        expected_lines = ["##MuTect=foo normal_sample_name=normal_sample tumor_sample_name=tumor_sample foo=bar", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL"]
        self.assertEquals(expected_lines, writer.lines())

    def test_normalize_doesnt_change_column_headers(self):
        writer = MockWriter()
        content = ["##MuTect=foo", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor_sample\tnormal_sample"]
        reader = MockReader(content)

        self.assertRaises(JQException, self.caller.normalize, writer, [reader])


#     #TODO: (cgates/kmeng) Remove when switched to new VcfRecord
#     def test_format_rounds(self):
#         tag = mutect.AlleleFreqTag()
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.2".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.2'), ('JQ_AF_MT', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.20".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.20'), ('JQ_AF_MT', '0.20')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204'), ('JQ_AF_MT', '0.2')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.205".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.205'), ('JQ_AF_MT', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.206".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.206'), ('JQ_AF_MT', '0.21')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:1.0".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '1.0'), ('JQ_AF_MT', '1.0')]), tag.format("alt", "filter", "", format_dict, 0))
#         
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:1.00".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '1.00'), ('JQ_AF_MT', '1.00')]), tag.format("alt", "filter", "", format_dict, 0))
# 
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204,0.3807".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204,0.3807'), ('JQ_AF_MT', '0.2,0.38')]), tag.format("alt", "filter", "", format_dict, 0))
#  
#         format_dict = OrderedDict(zip("A:FA".split(":"), "1:0.204,0.3807,0.2784".split(":")))
#         self.assertEqual(OrderedDict([('A', '1'), ('FA', '0.204,0.3807,0.2784'), ('JQ_AF_MT', '0.2,0.38,0.28')]), tag.format("alt", "filter", "", format_dict, 0))
