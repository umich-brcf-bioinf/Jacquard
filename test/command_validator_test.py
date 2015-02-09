#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=global-statement, unused-argument
#pylint: disable=too-many-instance-attributes
from __future__ import absolute_import
from argparse import Namespace
from testfixtures import TempDirectory
import jacquard.command_validator as command_validator
import jacquard.utils as utils
import os
import shutil
import subprocess
import sys
import test.mock_module as mock_module
import test.test_case as test_case
import time


GET_TEMP_WORKING_DIR_CALLED = False
CHECK_INPUT_EXISTS_CALLED = False
CHECK_INPUT_READABLE_CALLED = False
CHECK_INPUT_CORRECT_TYPE_CALLED = False
CHECK_OUTPUT_EXISTS_CALLED = False
CHECK_OUTPUT_CORRECT_TYPE_CALLED = False
CREATE_TEMP_WORKING_DIR_CALLED = False
CHECK_THERE_WILL_BE_OUTPUT_CALLED = False
CHECK_OVERWRITE_EXISTING_FILES_CALLED = False


class MakeDirTestCase(test_case.JacquardBaseTestCase):
    def test_makepath_singleLeaf(self):
        with TempDirectory() as output_dir:
            dest = os.path.join(output_dir.path, "foo_bar")
            command_validator._makepath(dest)
            self.assertTrue(os.path.isdir(dest))

    def test_makepath_okIfExists(self):
        with TempDirectory() as output_dir:
            dest = output_dir.path
            command_validator._makepath(dest)
            self.assertTrue(os.path.isdir(dest))

    def test_makepath_createsParents(self):
        with TempDirectory() as output_dir:
            dest = os.path.join(output_dir.path, "foo/bar")
            command_validator._makepath(dest)
            self.assertTrue(os.path.isdir(dest))

    def test_makepath_raisesOnWrongType(self):
        with TempDirectory() as output_dir:
            output_dir.write('foo', b'some text')
            dest = os.path.join(output_dir.path, "foo")
            self.assertRaises(OSError, command_validator._makepath, dest)

    def test_makepath_raisesIfExistingParentWrongType(self):
        with TempDirectory() as output_dir:
            output_dir.write('foo', b'some text')
            dest = os.path.join(output_dir.path, "foo/bar")
            self.assertRaises(OSError, command_validator._makepath, dest)

class CommandValidatorCheckOutputTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        self.original_check_output_correct_type = command_validator._check_output_correct_type
        command_validator._check_output_correct_type = mock_check_output_correct_type

    def tearDown(self):
        command_validator._check_output_correct_type = self.original_check_output_correct_type
        reset_called_methods()

    def test_check_output_exists_fileExists(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path, "output_file.vcf")

            command_validator._check_output_exists("awesomeCommand", output_file, "foo")

            self.assertTrue(CHECK_OUTPUT_CORRECT_TYPE_CALLED)

    def test_check_output_exists_fileDoesNotExist(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path, "output_file.vcf")

            command_validator._check_output_exists("awesomeCommand", output_file, "foo")
            self.assertTrue(True)

    def test_check_output_exists_directoryExists(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("output_dir")
            expected_output_dir = os.path.join(output_dir.path, "output_dir")

            command_validator._check_output_exists("awesomeCommand", expected_output_dir, "foo")
            self.assertTrue(CHECK_OUTPUT_CORRECT_TYPE_CALLED)

class CommandValidatorTestCase(test_case.JacquardBaseTestCase):
    def test_get_temp_working_dir_insideDestIfDestDirExists(self):
        with TempDirectory() as original_output_dir:
            original_output = os.path.join(original_output_dir.path, "dest_dir")
            os.mkdir(original_output)
            (temp_dir,
             output) = command_validator._get_temp_working_dir(original_output,
                                                               "directory")
        self.assertRegexpMatches(temp_dir,
                                 os.path.join(original_output,
                                              r"jacquard\.\d+\.\d+\.tmp"))
        self.assertEquals(temp_dir, output)

    def test_get_temp_working_dir_besideDestIfDestDirNotExists(self):
        with TempDirectory() as original_parent_dir:
            parent_dir = original_parent_dir.path
            original_output = os.path.join(parent_dir, "dest_dir")
            (temp_dir,
             output) = command_validator._get_temp_working_dir(original_output,
                                                               "directory")

        expected_temp_dir = os.path.join(parent_dir,
                                         r"jacquard\.\d+\.\d+\.tmp")
        self.assertRegexpMatches(temp_dir, expected_temp_dir)
        expected_output = os.path.join(parent_dir,
                                       r"jacquard\.\d+\.\d+\.tmp",
                                       "dest_dir")
        self.assertRegexpMatches(output, expected_output)

    def test_get_temp_working_dir_besideDestIfDestIsFile(self):
        with TempDirectory() as original_parent_dir:
            parent_dir = original_parent_dir.path
            original_output = os.path.join(parent_dir, "dest_file")
            (temp_dir,
             output) = command_validator._get_temp_working_dir(original_output,
                                                               "file")

        expected_temp_dir = os.path.join(parent_dir,
                                         r"jacquard\.\d+\.\d+\.tmp")
        self.assertRegexpMatches(temp_dir, expected_temp_dir)
        expected_output = os.path.join(parent_dir,
                                       r"jacquard\.\d+\.\d+\.tmp",
                                       "dest_file")
        self.assertRegexpMatches(output, expected_output)

    def test_get_temp_working_dir_distinct(self):
        orig_output_path = os.path.join("foo", "bar")
        temp_dir1, dummy = command_validator._get_temp_working_dir(orig_output_path,
                                                                   "file")
        time.sleep(0.01)
        temp_dir2, dummy = command_validator._get_temp_working_dir(orig_output_path,
                                                                   "file")
        self.assertNotEqual(temp_dir1, temp_dir2)

    def test_get_temp_working_dir_complexPathOk(self):
        original_output = os.path.join("foo.bar", "baz.hoopy", "frood")
        temp_dir, output = command_validator._get_temp_working_dir(original_output,
                                                                   "directory")
        self.assertRegexpMatches(temp_dir,
                                 os.path.join("foo.bar",
                                              "baz.hoopy",
                                              r"jacquard\.\d+\.\d+\.tmp"))
        self.assertRegexpMatches(output,
                                 os.path.join("foo.bar",
                                              "baz.hoopy",
                                              r"jacquard\.\d+\.\d+\.tmp",
                                              "frood"))

    def test_get_temp_working_dir_trailingSlashOk(self):
        original_output = os.path.join("foo", "bar") + os.sep
        temp_dir, output = command_validator._get_temp_working_dir(original_output,
                                                                   "directory")
        self.assertRegexpMatches(temp_dir,
                                 os.path.join("foo", "jacquard.*tmp"))
        self.assertRegexpMatches(output,
                                 os.path.join("foo", "jacquard.*tmp", "bar"))

    def test_get_temp_working_dir_currentDirOk(self):
        original_cwd = os.getcwd()
        try:
            with TempDirectory() as cwd:
                cwd_absolute_path = os.path.abspath(cwd.path)
                os.chdir(cwd_absolute_path)
                (temp_dir,
                 output) = command_validator._get_temp_working_dir(os.path.join("."),
                                                                   "directory")
                os.chdir(original_cwd)
                self.assertIn(cwd_absolute_path, temp_dir)
                self.assertEquals(temp_dir, output)
        finally:
            os.chdir(original_cwd)

    def test_get_temp_working_dir_relativePathMadeAbsolute(self):
        original_cwd = os.getcwd()
        try:
            with TempDirectory() as cwd:
                cwd_absolute_path = os.path.abspath(cwd.path)
                os.chdir(cwd_absolute_path)
                actual_dir, dummy = command_validator._get_temp_working_dir(os.path.join(".", "foo", "bar"),
                                                                            "directory")
                os.chdir(original_cwd)
                self.assertIn(cwd_absolute_path, actual_dir)
        finally:
            os.chdir(original_cwd)

    def test_check_input_exists_fileExists(self):
        with TempDirectory() as input_dir:
            input_dir.write("existant_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "existant_file.vcf")
            command_validator._check_input_exists(input_file)
            self.assertTrue(1 == 1)

    def test_check_input_exists_fileNonExists(self):
        with TempDirectory() as input_dir:
            input_file = os.path.join(input_dir.path, "nonexistant_file.vcf")
            self.assertRaisesRegexp(utils.UsageError,
                                    "Specified input .* does not exist",
                                    command_validator._check_input_exists, input_file)

    def test_check_input_exists_directoryExists(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("existant_dir")
            input_file = os.path.join(input_dir.path, "existant_dir")
            command_validator._check_input_exists(input_file)
            self.assertTrue(1 == 1)

    def test_check_input_exists_directoryNonExists(self):
        with TempDirectory() as input_dir:
            input_file = os.path.join(input_dir.path, "nonexistant_dir")
            self.assertRaisesRegexp(utils.UsageError,
                                    "Specified input .* does not exist",
                                    command_validator._check_input_exists,
                                    input_file)

    def test_check_input_readable_fileReadable(self):
        with TempDirectory() as input_dir:
            input_dir.write("readable_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "readable_file.vcf")
            command_validator._check_input_readable(input_file)

            self.assertTrue(1 == 1)

    def test_check_input_readable_fileNotReadable(self):
        with TempDirectory() as input_dir:
            input_dir.write("readable_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "readable_file.vcf")

            try:
                make_unreadable(input_file)
                self.assertRaisesRegexp(utils.UsageError,
                                        "Specified input .* cannot be read",
                                        command_validator._check_input_readable, input_file)

            finally:
                cleanup_unreadable(input_file)

    def test_check_input_readable_directoryReadable(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("readable_directory")
            readable_dir = os.path.join(input_dir.path, "readable_directory")
            command_validator._check_input_readable(readable_dir)

            self.assertTrue(1 == 1)

    def test_check_input_readable_directoryNotReadable(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("unreadable_directory")
            unreadable_dir = os.path.join(input_dir.path, "unreadable_directory")

            try:
                make_unreadable(unreadable_dir)
                self.assertRaisesRegexp(utils.UsageError,
                                        "Specified input .* cannot be read",
                                        command_validator._check_input_readable, unreadable_dir)

            finally:
                cleanup_unreadable(unreadable_dir)

    def test_check_input_correct_type_correct(self):
        with TempDirectory() as input_dir:
            correct_type = "directory"
            command_validator._check_input_correct_type("awesomeCommand", input_dir.path, correct_type)
            self.assertTrue(1 == 1)

            input_dir.write("input_file.vcf", "foo")
            correct_type = "file"
            input_file = os.path.join(input_dir.path, "input_file.vcf")
            command_validator._check_input_correct_type("awesomeCommand", input_file, correct_type)
            self.assertTrue(1 == 1)

    def test_check_input_correct_type_incorrect(self):
        with TempDirectory() as input_dir:
            input_dir.write("input_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "input_file.vcf")
            correct_type = "directory"
            self.assertRaisesRegexp(utils.UsageError,
                                    "The awesomeCommand command requires a directory.*input_file.vcf.*is a file",
                                    command_validator._check_input_correct_type,
                                    "awesomeCommand",
                                    input_file,
                                    correct_type)

            correct_type = "file"
            self.assertRaisesRegexp(utils.UsageError,
                                    "The awesomeCommand command requires a file.*" + input_dir.path + ".*is a directory",
                                    command_validator._check_input_correct_type,
                                    "awesomeCommand",
                                    input_dir.path,
                                    correct_type)

    def test_check_output_correct_type_correct(self):
        with TempDirectory() as output_dir:
            correct_type = "directory"
            command_validator._check_output_correct_type("awesomeCommand",
                                                         output_dir.path,
                                                         correct_type)
            self.assertTrue(1 == 1)

            output_dir.write("output_file.vcf", "foo")
            correct_type = "file"
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            command_validator._check_output_correct_type("awesomeCommand",
                                                         output_file,
                                                         correct_type)
            self.assertTrue(1 == 1)

    def test_check_output_correct_type_incorrect(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            correct_type = "directory"
            self.assertRaisesRegexp(utils.UsageError,
                                    "The awesomeCommand command outputs a directory.*output_file.vcf.*is a file",
                                    command_validator._check_output_correct_type,
                                    "awesomeCommand",
                                    output_file,
                                    correct_type)

            correct_type = "file"
            self.assertRaisesRegexp(utils.UsageError,
                                    "The awesomeCommand command outputs a file.*" + output_dir.path + ".*is a directory",
                                    command_validator._check_output_correct_type,
                                    "awesomeCommand",
                                    output_dir.path,
                                    correct_type)

    def test_create_temp_working_dir_noProblemIfTempAlreadyExists(self):
        with TempDirectory() as output_dir:
            tmp_dir = os.path.join(output_dir.path, "jacquard.tmp")
            os.mkdir(tmp_dir)
            args = Namespace(output=os.path.join(output_dir.path, "dest_dir"),
                             temp_working_dir=tmp_dir,
                             required_output_type="directory")

            command_validator._create_temp_working_dir(args)
            self.assertTrue(1 == 1)

    def test_create_temp_working_dir_createsIfNotExists(self):
        with TempDirectory() as output_dir:
            tmp_dir = os.path.join(output_dir.path, "jacquard.tmp")
            args = Namespace(output=os.path.join(output_dir.path, "dest_dir"),
                             temp_working_dir=tmp_dir,
                             required_output_type="directory")

            command_validator._create_temp_working_dir(args)
            self.assertTrue(os.path.isdir(tmp_dir))

    ## cgates: This test is a good idea, but doesn't work reliably on Windows 7.
    ## Can we remove it with a clean conscience? What's a better way?
    def xtest_create_temp_working_dir_doesNotExistCannotCreateTmp_directoryGivenAsOutput(self):
        with TempDirectory() as output_dir:
            unwriteable_dir = os.path.join(output_dir.path, "unwriteable_dir")
            output_file = os.path.join(unwriteable_dir, "output_file.vcf")

            try:
                make_unwriteable_dir(unwriteable_dir)
                self.assertRaisesRegexp(utils.UsageError,
                                        "Output directory .* cannot be read.",
                                        command_validator._create_temp_working_dir,
                                        output_file)
            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

    def test_create_temp_working_dir_doesExist_directoryGivenAsOutput(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("jacquard_tmp")
            args = Namespace(output=os.path.join(output_dir.path, "dest_dir"),
                             temp_working_dir=output_dir.path,
                             required_output_type="directory")

            command_validator._create_temp_working_dir(args)
            self.assertTrue(1 == 1)

    def test_check_overwrite_existing_files_raiseIfSingleFile(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            output_dir.write("output_file.vcf", "foo")
            self.assertRaisesRegexp(utils.UsageError,
                                    ("The foo command would "
                                     "overwrite the existing file "
                                     r"\[output_file.vcf\]; review "
                                     "command/output dir to avoid overwriting or "
                                     "use the flag '--force'."),
                                    command_validator._check_overwrite_existing_files,
                                    output_file,
                                    set([os.path.basename(output_file)]),
                                    "foo")

    def test_check_overwrite_existing_files_willNotRaiseIfNoOverwrite(self):
        with TempDirectory() as output_dir:
            output = os.path.join(output_dir.path, "foo.vcf")
            command_validator._check_overwrite_existing_files(output,
                                                              set(["foo.vcf"]),
                                                              "foo")
            self.assertTrue(1 == 1)

    def test_check_overwrite_existing_files_willOverwriteFiles(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file1.vcf", "foo")
            output_dir.write("output_file2.vcf", "foo")
            self.assertRaisesRegexp(utils.UsageError,
                                    "The foo command would overwrite 2 "
                                    r"existing files \[output_file1.vcf, output_file2.vcf\]; "
                                    "review command/output dir to avoid overwriting or "
                                    "use the flag '--force'.",
                                    command_validator._check_overwrite_existing_files,
                                    output_dir.path,
                                    set(["output_file1.vcf", "output_file2.vcf"]),
                                    "foo")


    def test_check_overwrite_existing_files_willOverwriteMoreThanFiveFiles(self):
        with TempDirectory() as output_dir:
            predicted_outputs = set()
            for x in xrange(1, 7):
                file_name = "out{}.vcf".format(str(x))
                output_dir.write(file_name, "foo")
                predicted_outputs.add(file_name)
            self.assertRaisesRegexp(utils.UsageError,
                                    "The foo command would overwrite 6 "
                                    r"existing files \[out1.vcf, out2.vcf, out3.vcf, out4.vcf, out5.vcf, "
                                    r"...\(1 file\(s\) omitted\)\]; "
                                    "review command/output dir to avoid "
                                    "overwriting or use the flag '--force'.",
                                    command_validator._check_overwrite_existing_files,
                                    output_dir.path,
                                    predicted_outputs,
                                    "foo")


    def test_check_overwrite_existing_files_willNotOverwriteDirectory(self):
        with TempDirectory() as output_dir:
            command_validator._check_overwrite_existing_files(output_dir.path,
                                                              set(["output_file1.vcf", "output_file2.vcf"]),
                                                              "foo")
            self.assertTrue(1 == 1)

    def test_check_overwrite_existing_files_forceFlagAdded(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file1.vcf", "foo")
            output_dir.write("output_file2.vcf", "foo")
            command_validator._check_overwrite_existing_files(output_dir.path,
                                                              set(["output_file1.vcf", "output_file2.vcf"]),
                                                              "foo",
                                                              1)
            self.assertTrue(True)

    def test_check_there_will_be_output(self):
        command_validator._check_there_will_be_output("awesomeCommand",
                                                      "input_dir",
                                                      ["out1", "out2"])
        self.assertTrue(True)

    def test_check_there_will_be_output_raisesIfNoPredictedOutput(self):
        predicted_outputs = []
        expected_regexp = ("Executing the awesomeCommand command with the "
                           r"input \[input_dir\] would not create any output "
                           "files. Review inputs and try again.")
        self.assertRaisesRegexp(utils.UsageError,
                                expected_regexp,
                                command_validator._check_there_will_be_output,
                                "awesomeCommand",
                                "input_dir",
                                predicted_outputs)


class CommandValidatorPreflightTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        self.original_get_temp_working_dir = command_validator._get_temp_working_dir
        self.original_check_input_exists = command_validator._check_input_exists
        self.original_check_input_readable = command_validator._check_input_readable
        self.original_check_input_correct_type = command_validator._check_input_correct_type
        self.original_check_output_exists = command_validator._check_output_exists
        self.original_check_output_correct_type = command_validator._check_output_correct_type
        self.original_create_temp_working_dir = command_validator._create_temp_working_dir
        self.original_check_there_will_be_output = command_validator._check_there_will_be_output
        self.original_check_overwrite_existing_files = command_validator._check_overwrite_existing_files

        command_validator._get_temp_working_dir = mock_get_temp_working_dir
        command_validator._check_input_exists = mock_check_input_exists
        command_validator._check_input_readable = mock_check_input_readable
        command_validator._check_input_correct_type = mock_check_input_correct_type
        command_validator._check_output_exists = mock_check_output_exists
        command_validator._check_output_correct_type = mock_check_output_correct_type
        command_validator._create_temp_working_dir = mock_create_temp_working_dir
        command_validator._check_there_will_be_output = mock_check_there_will_be_output
        command_validator._check_overwrite_existing_files = mock_check_overwrite_existing_files

    def tearDown(self):
        command_validator._get_temp_working_dir = self.original_get_temp_working_dir
        command_validator._check_input_exists = self.original_check_input_exists
        command_validator._check_input_readable = self.original_check_input_readable
        command_validator._check_input_correct_type = self.original_check_input_correct_type
        command_validator._check_output_exists = self.original_check_output_exists
        command_validator._check_output_correct_type = self.original_check_output_correct_type
        command_validator._create_temp_working_dir = self.original_create_temp_working_dir
        command_validator._check_there_will_be_output = self.original_check_there_will_be_output
        command_validator._check_overwrite_existing_files = self.original_check_overwrite_existing_files
        reset_called_methods()

    def test_preflight_inputFile(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_file = os.path.join(input_dir.path, "input.txt")
            output_file = os.path.join(output_dir.path, "output.txt")
            input_dir.write("input.txt", "foo")

            args = Namespace(input=input_file,
                             output=output_file,
                             subparser_name="foo",
                             force=0)

            command_validator.preflight(args, mock_module)

            self.assertTrue(CHECK_INPUT_EXISTS_CALLED)
            self.assertTrue(CHECK_INPUT_READABLE_CALLED)
            self.assertTrue(CHECK_INPUT_CORRECT_TYPE_CALLED)
            self.assertTrue(CHECK_OUTPUT_EXISTS_CALLED)
            self.assertTrue(not CHECK_OUTPUT_CORRECT_TYPE_CALLED)
            self.assertTrue(CREATE_TEMP_WORKING_DIR_CALLED)

    def test_preflight_outputFileExists(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_file = os.path.join(input_dir.path, "input.txt")
            input_dir.write("input.txt", "foo")
            output_dir.write("output.txt", "foo")
            output_file = os.path.join(output_dir.path, "output.txt")

            args = Namespace(input=input_file,
                             output=output_file,
                             subparser_name="foo",
                             force=0)

            command_validator.preflight(args, mock_module)

            self.assertTrue(GET_TEMP_WORKING_DIR_CALLED)
            self.assertTrue(CHECK_INPUT_EXISTS_CALLED)
            self.assertTrue(CHECK_INPUT_READABLE_CALLED)
            self.assertTrue(CHECK_INPUT_CORRECT_TYPE_CALLED)
            self.assertTrue(CHECK_OUTPUT_EXISTS_CALLED)
            self.assertTrue(CHECK_OUTPUT_CORRECT_TYPE_CALLED)
            self.assertTrue(CHECK_OVERWRITE_EXISTING_FILES_CALLED)

def mock_get_temp_working_dir(unused1, unused2):
    global GET_TEMP_WORKING_DIR_CALLED
    GET_TEMP_WORKING_DIR_CALLED = True
    return (None, None)

def mock_check_input_exists(unused):
    global CHECK_INPUT_EXISTS_CALLED
    CHECK_INPUT_EXISTS_CALLED = True

def mock_check_input_readable(unused):
    global CHECK_INPUT_READABLE_CALLED
    CHECK_INPUT_READABLE_CALLED = True

def mock_check_input_correct_type(unused1, unused2, unused3):
    global CHECK_INPUT_CORRECT_TYPE_CALLED
    CHECK_INPUT_CORRECT_TYPE_CALLED = True

def mock_check_output_exists(module_name, output, required_type):
    global CHECK_OUTPUT_EXISTS_CALLED
    CHECK_OUTPUT_EXISTS_CALLED = True
    if os.path.exists(output):
        mock_check_output_correct_type(module_name, output, required_type)

def mock_check_output_correct_type(unused1, unused2, unused3):
    global CHECK_OUTPUT_CORRECT_TYPE_CALLED
    CHECK_OUTPUT_CORRECT_TYPE_CALLED = True

def mock_create_temp_working_dir(unused):
    global CREATE_TEMP_WORKING_DIR_CALLED
    CREATE_TEMP_WORKING_DIR_CALLED = True

def mock_check_there_will_be_output(*args):
    global CHECK_THERE_WILL_BE_OUTPUT_CALLED
    CHECK_THERE_WILL_BE_OUTPUT_CALLED = True

def mock_check_overwrite_existing_files(unused1, unused2, unused3, unused4):
    global CHECK_OVERWRITE_EXISTING_FILES_CALLED
    CHECK_OVERWRITE_EXISTING_FILES_CALLED = True

def reset_called_methods():
    global GET_TEMP_WORKING_DIR_CALLED
    global CHECK_INPUT_EXISTS_CALLED
    global CHECK_INPUT_READABLE_CALLED
    global CHECK_INPUT_CORRECT_TYPE_CALLED
    global CHECK_OUTPUT_EXISTS_CALLED
    global CHECK_OUTPUT_CORRECT_TYPE_CALLED
    global CREATE_TEMP_WORKING_DIR_CALLED
    global CHECK_THERE_WILL_BE_OUTPUT_CALLED
    global CHECK_OVERWRITE_EXISTING_FILES_CALLED

    GET_TEMP_WORKING_DIR_CALLED = False
    CHECK_INPUT_EXISTS_CALLED = False
    CHECK_INPUT_READABLE_CALLED = False
    CHECK_INPUT_CORRECT_TYPE_CALLED = False
    CHECK_OUTPUT_EXISTS_CALLED = False
    CHECK_OUTPUT_CORRECT_TYPE_CALLED = False
    CREATE_TEMP_WORKING_DIR_CALLED = False
    CHECK_THERE_WILL_BE_OUTPUT_CALLED = False
    CHECK_OVERWRITE_EXISTING_FILES_CALLED = False

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

