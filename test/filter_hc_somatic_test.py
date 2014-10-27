import glob
import os
import unittest
import shutil
from StringIO import StringIO
import sys
import testfixtures
from testfixtures import TempDirectory
from jacquard.filter_hc_somatic import filter_somatic_positions, write_somatic, find_somatic_positions

import jacquard.filter_hc_somatic as filter_hc_somatic
import jacquard.logger as logger

VCF_HEADER="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB\n"

mock_log_called = False
  
def mock_log(msg, *args):
    global mock_log_called
    mock_log_called = True
    
class FilterSomaticTestCase(unittest.TestCase):
    def setUp(self):
        self.output = StringIO()
        self.saved_stderr = sys.stderr
        sys.stderr = self.output
        self.original_info = logger.info
        self.original_error = logger.error
        self.original_warning = logger.warning
        self.original_debug = logger.debug
        self._change_mock_logger()
        
    def tearDown(self):
        self.output.close()
        sys.stderr = self.saved_stderr
        self._reset_mock_logger()
    
    def _change_mock_logger(self):
        global mock_log_called
        mock_log_called = False
        global mock_log
        logger.info = mock_log
        logger.error = mock_log
        logger.warning = mock_log
        logger.debug = mock_log
        
    def _reset_mock_logger(self):
        logger.info = self.original_info
        logger.error = self.original_error
        logger.warning = self.original_warning
        logger.debug = self.original_debug
        
    def test_findSomaticPositions(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.snp.vcf","##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP:JQ_HC_SOM_VS\t234:1\n")
            input_dir.write("A.indel.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2353\t.\tA\tGT\t.\t.\tfoo\tDP:JQ_HC_SOM_VS\t234:1\n")
            
            file1 = os.path.join(input_dir.path, "A.snp.vcf")
            file2 = os.path.join(input_dir.path, "A.indel.vcf")
            
            somatic_positions, somatic_positions_header  = find_somatic_positions([file1, file2], output_dir.path)
            self.assertEqual({'1^2353': 1, '1^2352': 1}, somatic_positions)
            self.assertEqual("##jacquard.filterHCSomatic.total_highConfidence_somatic_positions=2\n", somatic_positions_header)
            
            input_dir.cleanup()
            output_dir.cleanup()
    
    def test_findSomaticPositions_invalidInput(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.snp.vcf","##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234:1\n")
            input_dir.write("A.indel.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2353\t.\tA\tGT\t.\t.\tfoo\tDP:JQ_SOM_VS\t234:1\n")
            
            file1 = os.path.join(input_dir.path, "A.snp.vcf")
            file2 = os.path.join(input_dir.path, "A.indel.vcf")
            
            with self.assertRaises(SystemExit) as cm:
                somatic_positions = find_somatic_positions([file1, file2], output_dir.path)
            self.assertEqual(cm.exception.code, 1)
            
            input_dir.cleanup()
            output_dir.cleanup()
            
    def test_writeSomatic_outputFileForEachInputFile(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input1 = input_dir.write('mutect.vcf', "##MuTect\n"+VCF_HEADER)
            input2 = input_dir.write('varscan.vcf', "##source=VarScan2\n"+VCF_HEADER)
            
            in_files = [input1, input2]
            somatic_positions = {}
            execution_context = ["##foo", "##bar"]
            excluded_variants = write_somatic(in_files, output_dir.path, somatic_positions, execution_context)
            self.assertEqual(["mutect_HCsomatic.vcf", "varscan_HCsomatic.vcf"], output_dir.actual())

    
    def Ytest_writeSomatic5(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/reference_files/filter_somatic_test/input/"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
        output_dir = script_dir + "/reference_files/filter_somatic_test/output"
        try:
            os.mkdir(output_dir)
        except:
            pass
        somatic_positions = {'1^2353': 1, '1^2352': 1}
        execution_context = ["##foo.version=0.1", "##bar"]
        excluded_variants = write_somatic(in_files, output_dir, somatic_positions, execution_context)
        
        self.assertEqual("##jacquard.filterHCSomatic.excluded_variants=5\n", excluded_variants)
        shutil.rmtree(output_dir)
        
    def Ytest_writeSomatic32(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/reference_files/filter_somatic_test/input/"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
        output_dir = script_dir + "/reference_files/filter_somatic_test/output"
        try:
            os.mkdir(output_dir)
        except:
            pass
        somatic_positions = {'1^15996': 1, '1^2352': 1}
        execution_context = ["##foo.version=0.1", "##bar"]
        excluded_variants = write_somatic(in_files, output_dir, somatic_positions, execution_context)
        
        self.assertEqual("##jacquard.filterHCSomatic.excluded_variants=32\n", excluded_variants)
        shutil.rmtree(output_dir)

            
