import unittest
from jacquard.vcf_record import VcfRecord

class VcfRecordTestCase(unittest.TestCase):

    def testInit(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|SAMPLE1|SAMPLE2\n".replace('|',"\t")
        record = VcfRecord(input_line)
        self.assertEquals("CHROM", record.chrom)
        self.assertEquals("POS", record.pos)
        self.assertEquals("ID", record.id)
        self.assertEquals("REF", record.ref)
        self.assertEquals("ALT", record.alt)
        self.assertEquals("QUAL", record.qual)
        self.assertEquals("FILTER", record.filter)
        self.assertEquals("INFO", record.info)
        self.assertEquals("FORMAT", record.format)
        self.assertEquals(["SAMPLE1","SAMPLE2"], record.samples)
        
    def testInit_format_set(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord(input_line)
        self.assertEquals(set(["F1","F2","F3"]), record.format_set)

    def testInit_sample_dict(self):
        input_line = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|F1:F2:F3|SA.1:SA.2:SA.3|SB.1:SB.2:SB.3\n".replace('|',"\t")
        record = VcfRecord(input_line)
        self.assertEquals([0, 1], record.sample_dict.keys())
        self.assertEquals({"F1":"SA.1","F2":"SA.2","F3":"SA.3"}, record.sample_dict[0])
        self.assertEquals({"F1":"SB.1","F2":"SB.2","F3":"SB.3"}, record.sample_dict[1])
        
