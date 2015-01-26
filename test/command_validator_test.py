#pylint: disable=line-too-long, too-many-public-methods, invalid-name
from __future__ import absolute_import

import os
import shutil
import subprocess
import sys

from testfixtures import TempDirectory

import jacquard.command_validator as command_validator
import jacquard.utils as utils
import test.test_case as test_case

class CommandValidatorTestCase(test_case.JacquardBaseTestCase):
    def test_check_input_exists_fileExistant(self):
        with TempDirectory() as input_dir:
            input_dir.write("existant_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "existant_file.vcf")
            command_validator.check_input_exists(input_file)
            self.assertTrue(1==1)

    def test_check_input_exists_fileNonExistant(self):
        with TempDirectory() as input_dir:
            input_file = os.path.join(input_dir.path, "nonexistant_file.vcf")
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified input .* does not exist",
                                    command_validator.check_input_exists, input_file)

    def test_check_input_exists_directoryExistant(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("existant_dir")
            input_file = os.path.join(input_dir.path, "existant_dir")
            command_validator.check_input_exists(input_file)
            self.assertTrue(1==1)

    def test_check_input_exists_directoryNonExistant(self):
        with TempDirectory() as input_dir:
            input_file = os.path.join(input_dir.path, "nonexistant_dir")
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified input .* does not exist",
                                    command_validator.check_input_exists, input_file)

    def test_check_input_readable_fileReadable(self):
        with TempDirectory() as input_dir:
            input_dir.write("readable_file.vcf","foo")
            input_file = os.path.join(input_dir.path, "readable_file.vcf")
            command_validator.check_input_readable(input_file)

            self.assertTrue(1==1)

    def test_check_input_readable_fileNotReadable(self):
        with TempDirectory() as input_dir:
            input_dir.write("readable_file.vcf","foo")
            input_file = os.path.join(input_dir.path, "readable_file.vcf")

            try:
                make_unreadable(input_file)
                self.assertRaisesRegexp(utils.JQException,
                                        "Specified input .* cannot be read",
                                        command_validator.check_input_readable, input_file)

            finally:
                cleanup_unreadable(input_file)

    def test_check_input_readable_directoryReadable(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("readable_directory")
            readable_dir = os.path.join(input_dir.path, "readable_directory")
            command_validator.check_input_readable(readable_dir)

            self.assertTrue(1==1)

    def test_check_input_readable_directoryNotReadable(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("unreadable_directory")
            unreadable_dir = os.path.join(input_dir.path, "unreadable_directory")

            try:
                make_unreadable(unreadable_dir)
                self.assertRaisesRegexp(utils.JQException,
                                        "Specified input .* cannot be read",
                                        command_validator.check_input_readable, unreadable_dir)

            finally:
                cleanup_unreadable(unreadable_dir)

    def test_check_input_correct_type_correct(self):
        with TempDirectory() as input_dir:
            correct_type = "directory"
            command_validator.check_input_correct_type(input_dir.path, correct_type)
            self.assertTrue(1==1)

            input_dir.write("input_file.vcf","foo")
            correct_type = "file"
            input_file = os.path.join(input_dir.path,"input_file.vcf")
            command_validator.check_input_correct_type(input_file, correct_type)
            self.assertTrue(1==1)

    def test_check_input_correct_type_incorrect(self):
        with TempDirectory() as input_dir:
            input_dir.write("input_file.vcf","foo")
            input_file = os.path.join(input_dir.path,"input_file.vcf")
            correct_type = "directory"
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified input .* does not match command\'s required file type",
                                    command_validator.check_input_correct_type,
                                    input_file, correct_type)

            correct_type = "file"
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified input .* does not match command\'s required file type",
                                    command_validator.check_input_correct_type,
                                    input_dir.path, correct_type)

    def test_check_output_exists_fileExistant(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path,"output_file.vcf")

            command_validator.check_output_exists(output_file)
            self.assertTrue(1==1)

    def test_check_output_exists_fileNonExistantCanCreateOutput(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path,"output_file.vcf")

            command_validator.check_output_exists(output_file)
            self.assertTrue(1==1)

    def test_check_output_exists_fileNonExistantCannotCreateOutput(self):
        with TempDirectory() as output_dir:
            unwriteable_dir = os.path.join(output_dir.path, "unwriteable_dir")
            output_file = os.path.join(unwriteable_dir, "output_file.vcf")

            try:
                make_unwriteable_dir(unwriteable_dir)
                self.assertRaisesRegexp(utils.JQException,
                                        "Specified output .* does not exist and cannot be created",
                                        command_validator.check_output_exists,
                                        output_file)
            finally:
#                 pass
                cleanup_unwriteable_dir(unwriteable_dir)

    def test_check_output_exists_directoryExistant(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("output_dir")
            expected_output_dir = os.path.join(output_dir.path,"output_dir")

            command_validator.check_output_exists(expected_output_dir)
            self.assertTrue(1==1)

    def test_check_output_exists_directoryNonExistantCanCreateOutput(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path,"output_dir")

            command_validator.check_output_exists(output_file)
            self.assertTrue(1==1)

    def test_check_output_exists_directoryNonExistantCannotCreateOutput(self):
        with TempDirectory() as output_dir:
            unwriteable_dir = os.path.join(output_dir.path, "unwriteable_dir")
            output_file = os.path.join(unwriteable_dir,"output_dir")

            try:
                make_unwriteable_dir(unwriteable_dir)
                self.assertRaisesRegexp(utils.JQException,
                                        "Specified output .* does not exist and cannot be created",
                                        command_validator.check_output_exists,
                                        output_file)
            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

    def test_check_output_correct_type_correct(self):
        with TempDirectory() as output_dir:
            correct_type = "directory"
            command_validator.check_output_correct_type(output_dir.path, correct_type)
            self.assertTrue(1==1)

            output_dir.write("output_file.vcf","foo")
            correct_type = "file"
            output_file = os.path.join(output_dir.path,"output_file.vcf")
            command_validator.check_output_correct_type(output_file, correct_type)
            self.assertTrue(1==1)

    def test_check_output_correct_type_incorrect(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf","foo")
            output_file = os.path.join(output_dir.path,"output_file.vcf")
            correct_type = "directory"
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified output .* does not match command\'s required file type",
                                    command_validator.check_output_correct_type,
                                    output_file, correct_type)

            correct_type = "file"
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified output .* does not match command\'s required file type",
                                    command_validator.check_output_correct_type,
                                    output_dir.path, correct_type)

    def test_check_tmpdir_exists_doesExist_fileGivenAsOutput(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("jacquard_tmp")
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            command_validator.check_tmpdir_exists(output_file)
            self.assertTrue(1==1)

    def test_check_tmpdir_exists_doesNotExistCanCreateTmp_fileGivenAsOutput(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            command_validator.check_tmpdir_exists(output_file)
            self.assertTrue(1==1)

#TODO: (jebene) - fix this test!
    def xtest_check_tmpdir_exists_doesNotExistCannotCreateTmp_directoryGivenAsOutput(self):
        #pylint: disable=anomalous-backslash-in-string
        with TempDirectory() as output_dir:
            unwriteable_dir = os.path.join(output_dir.path, "unwriteable_dir")
            output_file = os.path.join(unwriteable_dir, "output_file.vcf")

            try:
                make_unwriteable_dir(unwriteable_dir)
                self.assertRaisesRegexp(utils.JQException,
                                        "Output directory .* cannot be read.",
                                         command_validator.check_tmpdir_exists, output_file)
            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

    def test_check_tmpdir_exists_doesExist_directoryGivenAsOutput(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("jacquard_tmp")
            command_validator.check_tmpdir_exists(output_dir.path)
            self.assertTrue(1==1)

    def test_check_overwrite_existing_files_willOverwriteFile(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            output_dir.write("output_file.vcf","foo")
            self.assertRaisesRegexp(utils.JQException,
                                    "ERROR: The command .* would "
                                    "overwrite existing files .*; review "
                                    "command/output dir to avoid overwriting or "
                                    "use the flag '--force'. Type 'jacquard -h' "
                                    "for more details",
                                    command_validator.check_overwrite_existing_files,
                                    output_file, set([os.path.basename(output_file)]),
                                    "foo")

    def test_check_overwrite_existing_files_willNotOverwriteFile(self):
        command_validator.check_overwrite_existing_files("foo.vcf",
                            set(["output_file1.vcf"]),
                            "foo")
        self.assertTrue(1==1)

    def test_check_overwrite_existing_files_willOverwriteDirectory(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file1.vcf","foo")
            output_dir.write("output_file2.vcf","foo")
            self.assertRaisesRegexp(utils.JQException,
                                    "ERROR: The command .* would "
                                    "overwrite existing files .*; review "
                                    "command/output dir to avoid overwriting or "
                                    "use the flag '--force'. Type 'jacquard -h' "
                                    "for more details",
                                    command_validator.check_overwrite_existing_files,
                                    output_dir.path,
                                    set(["output_file1.vcf","output_file2.vcf"]),
                                    "foo")

    def test_check_overwrite_existing_files_willNotOverwriteDirectory(self):
        with TempDirectory() as output_dir:
            command_validator.check_overwrite_existing_files(output_dir.path,
                                set(["output_file1.vcf","output_file2.vcf"]),
                                "foo")
            self.assertTrue(1==1)

    def test_check_overwrite_existing_files_forceFlagAdded(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file1.vcf","foo")
            output_dir.write("output_file2.vcf","foo")
            command_validator.check_overwrite_existing_files(output_dir.path,
                                                             set(["output_file1.vcf","output_file2.vcf"]),
                                                             "foo",
                                                             1)
            self.assertTrue(1==1)


def is_windows_os():
    return sys.platform.lower().startswith("win")

def make_unwriteable_dir(unwriteable_dir):
    os.mkdir(unwriteable_dir, 0555)

    if is_windows_os():
        FNULL = open(os.devnull, 'w')
        subprocess.call("icacls {0} /deny Everyone:W".format(unwriteable_dir), stdout=FNULL, stderr=subprocess.STDOUT)

def cleanup_unwriteable_dir(unwriteable_dir):
    if is_windows_os():
        FNULL = open(os.devnull, 'w')
        subprocess.call("icacls {0} /reset /t /c".format(unwriteable_dir), stdout=FNULL, stderr=subprocess.STDOUT)

    shutil.rmtree(unwriteable_dir)

def make_unreadable(unreadable_dir):
    os.chmod(unreadable_dir, 0333)

    if is_windows_os():
        FNULL = open(os.devnull, 'r')
        subprocess.call("icacls {0} /deny Everyone:R".format(unreadable_dir), stdout=FNULL, stderr=subprocess.STDOUT)

def cleanup_unreadable(unreadable):
    if is_windows_os():
        FNULL = open(os.devnull, 'r')
        subprocess.call("icacls {0} /reset /t /c".format(unreadable), stdout=FNULL, stderr=subprocess.STDOUT)

    if os.path.isdir(unreadable):
        os.rmdir(unreadable)
    elif os.path.isfile(unreadable):
        os.remove(unreadable)

