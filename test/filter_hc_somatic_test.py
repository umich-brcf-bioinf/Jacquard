#!/usr/bin/python2.7
import glob
import os
import unittest
import shutil
import testfixtures
from testfixtures import TempDirectory
from bin.filter_hc_somatic import filter_somatic_positions, write_somatic, find_somatic_positions

class FilterSomaticTestCase(unittest.TestCase):
    def test_findSomaticPositions(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_dir.write("A.snp.vcf","##source=VarScan2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP\t234\n1\t2352\t.\tA\tG\t.\t.\tfoo\tDP:JQ_HC_SOM_VS\t234:1\n")
            input_dir.write("A.indel.vcf","##source=VarScan2\n#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tSAMPLE\n1\t2353\t.\tA\tGT\t.\t.\tfoo\tDP:JQ_HC_SOM_VS\t234:1\n")
            
            file1 = os.path.join(input_dir.path, "A.snp.vcf")
            file2 = os.path.join(input_dir.path, "A.indel.vcf")
            
            somatic_positions = find_somatic_positions([file1, file2], output_dir.path)
            self.assertEqual({'1^2353': 1, '1^2352': 1}, somatic_positions)
            
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
            
    def test_writeSomatic(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/filter_somatic_test/input/"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
        output_dir = script_dir + "/filter_somatic_test/output"
        try:
            os.mkdir(output_dir)
        except:
            pass
        somatic_positions = {'1^2353': 1, '1^2352': 1}
        
        excluded_variants = write_somatic(in_files, output_dir, somatic_positions)
        
        self.assertEqual("##jacquard.filterHCSomatic.excluded_variants=5\n", excluded_variants)
        shutil.rmtree(output_dir)
        
    def test_writeSomatic(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_dir = script_dir + "/filter_somatic_test/input/"
        in_files = sorted(glob.glob(os.path.join(input_dir,"*.vcf")))
        output_dir = script_dir + "/filter_somatic_test/output"
        try:
            os.mkdir(output_dir)
        except:
            pass
        somatic_positions = {'1^15996': 1, '1^2352': 1}
        
        excluded_variants = write_somatic(in_files, output_dir, somatic_positions)
        
        self.assertEqual("##jacquard.filterHCSomatic.excluded_variants=14\n", excluded_variants)
        shutil.rmtree(output_dir)

            