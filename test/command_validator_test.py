#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=global-statement, unused-argument, blacklisted-name
#pylint: disable=too-many-instance-attributes
from __future__ import absolute_import

from argparse import Namespace
import os
import shutil
import subprocess
import sys

from testfixtures import TempDirectory

import jacquard.command_validator as command_validator
import jacquard.utils as utils
import test.mock_module as mock_module
import test.test_case as test_case

CHECK_INPUT_EXISTS_CALLED = False
CHECK_INPUT_READABLE_CALLED = False
CHECK_INPUT_CORRECT_TYPE_CALLED = False
CHECK_OUTPUT_EXISTS_CALLED = False
CHECK_CAN_CREATE_OUTPUT_CALLED = False
CHECK_OUTPUT_CORRECT_TYPE_CALLED = False
CHECK_TMPDIR_EXISTS_CALLED = False
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
        self.original_check_can_create_output = command_validator._check_can_create_output
        self.original_check_output_correct_type = command_validator._check_output_correct_type
        command_validator._check_can_create_output = mock_check_can_create_output
        command_validator._check_output_correct_type = mock_check_output_correct_type

    def tearDown(self):
        command_validator._check_can_create_output = self.original_check_can_create_output
        command_validator._check_output_correct_type = self.original_check_output_correct_type
        reset_called_methods()

    def test_check_output_exists_fileExistant(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path,"output_file.vcf")

            command_validator._check_output_exists(output_file, "foo")

            self.assertTrue(CHECK_OUTPUT_CORRECT_TYPE_CALLED)

    def test_check_output_exists_fileNonExistant(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path,"output_file.vcf")

            command_validator._check_output_exists(output_file, "foo")
            self.assertTrue(CHECK_CAN_CREATE_OUTPUT_CALLED)

    def test_check_output_exists_directoryExistant(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("output_dir")
            expected_output_dir = os.path.join(output_dir.path,"output_dir")

            command_validator._check_output_exists(expected_output_dir, "foo")
            self.assertTrue(CHECK_OUTPUT_CORRECT_TYPE_CALLED)

class CommandValidatorTestCase(test_case.JacquardBaseTestCase):
    def test_check_input_exists_fileExistant(self):
        with TempDirectory() as input_dir:
            input_dir.write("existant_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "existant_file.vcf")
            command_validator._check_input_exists(input_file)
            self.assertTrue(1==1)

    def test_check_input_exists_fileNonExistant(self):
        with TempDirectory() as input_dir:
            input_file = os.path.join(input_dir.path, "nonexistant_file.vcf")
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified input .* does not exist",
                                    command_validator._check_input_exists, input_file)

    def test_check_input_exists_directoryExistant(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("existant_dir")
            input_file = os.path.join(input_dir.path, "existant_dir")
            command_validator._check_input_exists(input_file)
            self.assertTrue(1==1)

    def test_check_input_exists_directoryNonExistant(self):
        with TempDirectory() as input_dir:
            input_file = os.path.join(input_dir.path, "nonexistant_dir")
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified input .* does not exist",
                                    command_validator._check_input_exists, input_file)

    def test_check_input_readable_fileReadable(self):
        with TempDirectory() as input_dir:
            input_dir.write("readable_file.vcf","foo")
            input_file = os.path.join(input_dir.path, "readable_file.vcf")
            command_validator._check_input_readable(input_file)

            self.assertTrue(1==1)

    def test_check_input_readable_fileNotReadable(self):
        with TempDirectory() as input_dir:
            input_dir.write("readable_file.vcf","foo")
            input_file = os.path.join(input_dir.path, "readable_file.vcf")

            try:
                make_unreadable(input_file)
                self.assertRaisesRegexp(utils.JQException,
                                        "Specified input .* cannot be read",
                                        command_validator._check_input_readable, input_file)

            finally:
                cleanup_unreadable(input_file)

    def test_check_input_readable_directoryReadable(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("readable_directory")
            readable_dir = os.path.join(input_dir.path, "readable_directory")
            command_validator._check_input_readable(readable_dir)

            self.assertTrue(1==1)

    def test_check_input_readable_directoryNotReadable(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("unreadable_directory")
            unreadable_dir = os.path.join(input_dir.path, "unreadable_directory")

            try:
                make_unreadable(unreadable_dir)
                self.assertRaisesRegexp(utils.JQException,
                                        "Specified input .* cannot be read",
                                        command_validator._check_input_readable, unreadable_dir)

            finally:
                cleanup_unreadable(unreadable_dir)

    def test_check_input_correct_type_correct(self):
        with TempDirectory() as input_dir:
            correct_type = "directory"
            command_validator._check_input_correct_type(input_dir.path, correct_type)
            self.assertTrue(1==1)

            input_dir.write("input_file.vcf","foo")
            correct_type = "file"
            input_file = os.path.join(input_dir.path,"input_file.vcf")
            command_validator._check_input_correct_type(input_file, correct_type)
            self.assertTrue(1==1)

    def test_check_input_correct_type_incorrect(self):
        with TempDirectory() as input_dir:
            input_dir.write("input_file.vcf","foo")
            input_file = os.path.join(input_dir.path,"input_file.vcf")
            correct_type = "directory"
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified input .* does not match command\'s required file type",
                                    command_validator._check_input_correct_type,
                                    input_file, correct_type)

            correct_type = "file"
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified input .* does not match command\'s required file type",
                                    command_validator._check_input_correct_type,
                                    input_dir.path, correct_type)

    def test_check_can_create_output_canCreateOutputParentAndFile(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path,"parent_dir/output_file.txt")

            command_validator._check_can_create_output(output_file)
            self.assertIn("parent_dir", os.listdir(output_dir.path))

    def test_check_can_create_output_relativeDirectory(self):
        current_working_dir = os.getcwd()
        with TempDirectory() as output_dir:
            try:
                os.chdir(output_dir.path)
                command_validator._check_can_create_output("target_dir")
                self.assertEquals([], os.listdir(output_dir.path))
            finally: 
                os.chdir(current_working_dir)

    def test_check_can_create_output_canCreateParentOutputDirectory(self):
        with TempDirectory() as output_dir:
            target_dir = os.path.join(output_dir.path,"nested_dir/target_dir")

            command_validator._check_can_create_output(target_dir)
            self.assertIn("nested_dir", os.listdir(output_dir.path))

    def test_check_can_create_output_cannotCreateOutputFile(self):
        with TempDirectory() as output_dir:
            unwriteable_dir = os.path.join(output_dir.path, "unwriteable_dir")
            output_file = os.path.join(unwriteable_dir, "nested_dir/output_file.vcf")

            try:
                make_unwriteable_dir(unwriteable_dir)
                self.assertRaisesRegexp(utils.JQException,
                                        "Specified output .* does not exist and cannot be created",
                                        command_validator._check_can_create_output,
                                        output_file)
            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

    def test_check_can_create_output_cannotCreateOutputDirectory(self):
        with TempDirectory() as output_dir:
            unwriteable_dir = os.path.join(output_dir.path, "unwriteable_dir")
            target_dir = os.path.join(unwriteable_dir,"nested_dir/target_dir")

            try:
                make_unwriteable_dir(unwriteable_dir)
                self.assertRaisesRegexp(utils.JQException,
                                        "Specified output .* does not exist and cannot be created",
                                        command_validator._check_can_create_output,
                                        target_dir)
            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

    def test_check_output_correct_type_correct(self):
        with TempDirectory() as output_dir:
            correct_type = "directory"
            command_validator._check_output_correct_type(output_dir.path, correct_type)
            self.assertTrue(1==1)

            output_dir.write("output_file.vcf","foo")
            correct_type = "file"
            output_file = os.path.join(output_dir.path,"output_file.vcf")
            command_validator._check_output_correct_type(output_file, correct_type)
            self.assertTrue(1==1)

    def test_check_output_correct_type_incorrect(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf","foo")
            output_file = os.path.join(output_dir.path,"output_file.vcf")
            correct_type = "directory"
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified output .* does not match command\'s required file type",
                                    command_validator._check_output_correct_type,
                                    output_file, correct_type)

            correct_type = "file"
            self.assertRaisesRegexp(utils.JQException,
                                    "Specified output .* does not match command\'s required file type",
                                    command_validator._check_output_correct_type,
                                    output_dir.path, correct_type)

    def test_check_tmpdir_exists_doesExist_fileGivenAsOutput(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("jacquard_tmp")
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            command_validator._check_tmpdir_exists(output_file)
            self.assertTrue(1==1)

    def test_check_tmpdir_exists_doesNotExistCanCreateTmp_fileGivenAsOutput(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            command_validator._check_tmpdir_exists(output_file)
            self.assertTrue(1==1)

    def xtest_check_tmpdir_exists_doesNotExistCannotCreateTmp_directoryGivenAsOutput(self):
        with TempDirectory() as output_dir:
            unwriteable_dir = os.path.join(output_dir.path, "unwriteable_dir")
            output_file = os.path.join(unwriteable_dir, "output_file.vcf")

            try:
                make_unwriteable_dir(unwriteable_dir)
                self.assertRaisesRegexp(utils.JQException,
                                        "Output directory .* cannot be read.",
                                         command_validator._check_tmpdir_exists, output_file)
            finally:
                cleanup_unwriteable_dir(unwriteable_dir)

    def test_check_tmpdir_exists_doesExist_directoryGivenAsOutput(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("jacquard_tmp")
            command_validator._check_tmpdir_exists(output_dir.path)
            self.assertTrue(1==1)

    def test_check_overwrite_existing_files_raiseIfFileWouldBeOverwritten(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            output_dir.write("output_file.vcf","foo")
            self.assertRaisesRegexp(utils.JQException,
                                    "ERROR: The command .* would "
                                    "overwrite existing files .*; review "
                                    "command/output dir to avoid overwriting or "
                                    "use the flag '--force'. Type 'jacquard -h' "
                                    "for more details",
                                    command_validator._check_overwrite_existing_files,
                                    output_file, set([os.path.basename(output_file)]),
                                    "foo")

    def test_check_overwrite_existing_files_willNotRaiseIfNoOverwrite(self):
        with TempDirectory() as output_dir:
            output = os.path.join(output_dir.path, "foo.vcf")
            command_validator._check_overwrite_existing_files(output,
                                                              set(["foo.vcf"]),
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
                                    command_validator._check_overwrite_existing_files,
                                    output_dir.path,
                                    set(["output_file1.vcf","output_file2.vcf"]),
                                    "foo")

    def test_check_overwrite_existing_files_willNotOverwriteDirectory(self):
        with TempDirectory() as output_dir:
            command_validator._check_overwrite_existing_files(output_dir.path,
                                set(["output_file1.vcf","output_file2.vcf"]),
                                "foo")
            self.assertTrue(1==1)

    def test_check_overwrite_existing_files_forceFlagAdded(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file1.vcf","foo")
            output_dir.write("output_file2.vcf","foo")
            command_validator._check_overwrite_existing_files(output_dir.path,
                                                             set(["output_file1.vcf","output_file2.vcf"]),
                                                             "foo",
                                                             1)
            self.assertTrue(1==1)

class CommandValidatorPreflightTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        self.original_check_input_exists = command_validator._check_input_exists
        self.original_check_input_readable = command_validator._check_input_readable
        self.original_check_input_correct_type = command_validator._check_input_correct_type
        self.original_check_output_exists = command_validator._check_output_exists
        self.original_check_can_create_output = command_validator._check_can_create_output
        self.original_check_output_correct_type = command_validator._check_output_correct_type
        self.original_check_tmpdir_exists = command_validator._check_tmpdir_exists
        self.original_check_overwrite_existing_files = command_validator._check_overwrite_existing_files

        command_validator._check_input_exists = mock_check_input_exists
        command_validator._check_input_readable = mock_check_input_readable
        command_validator._check_input_correct_type = mock_check_input_correct_type
        command_validator._check_output_exists = mock_check_output_exists
        command_validator._check_can_create_output = mock_check_can_create_output
        command_validator._check_output_correct_type = mock_check_output_correct_type
        command_validator._check_tmpdir_exists = mock_check_tmpdir_exists
        command_validator._check_overwrite_existing_files = mock_check_overwrite_existing_files

    def tearDown(self):
        command_validator._check_input_exists = self.original_check_input_exists
        command_validator._check_input_readable = self.original_check_input_readable
        command_validator._check_input_correct_type = self.original_check_input_correct_type
        command_validator._check_output_exists = self.original_check_output_exists
        command_validator._check_can_create_output = self.original_check_can_create_output
        command_validator._check_output_correct_type = self.original_check_output_correct_type
        command_validator._check_tmpdir_exists = self.original_check_tmpdir_exists
        command_validator._check_overwrite_existing_files = self.original_check_overwrite_existing_files
        reset_called_methods()

    def test_preflight_inputFile(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_file = os.path.join(input_dir.path, "input.txt")
            output_file = os.path.join(output_dir.path, "output.txt")
            input_dir.write("input.txt","foo")

            args = Namespace(input=input_file,
                             output=output_file,
                             subparser_name="foo",
                             force=0)

            command_validator.preflight(args, mock_module)

            self.assertTrue(CHECK_INPUT_EXISTS_CALLED)
            self.assertTrue(CHECK_INPUT_READABLE_CALLED)
            self.assertTrue(CHECK_INPUT_CORRECT_TYPE_CALLED)
            self.assertTrue(CHECK_OUTPUT_EXISTS_CALLED)
            self.assertTrue(CHECK_CAN_CREATE_OUTPUT_CALLED)
            self.assertTrue(not CHECK_OUTPUT_CORRECT_TYPE_CALLED)
            self.assertTrue(CHECK_TMPDIR_EXISTS_CALLED)

    def test_preflight_outputFileExists(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            input_file = os.path.join(input_dir.path, "input.txt")
            input_dir.write("input.txt","foo")
            output_dir.write("output.txt", "foo")
            output_file = os.path.join(output_dir.path, "output.txt")

            args = Namespace(input=input_file,
                             output=output_file,
                             subparser_name="foo",
                             force=0)

            command_validator.preflight(args, mock_module)

            self.assertTrue(CHECK_INPUT_EXISTS_CALLED)
            self.assertTrue(CHECK_INPUT_READABLE_CALLED)
            self.assertTrue(CHECK_INPUT_CORRECT_TYPE_CALLED)
            self.assertTrue(CHECK_OUTPUT_EXISTS_CALLED)
            self.assertTrue(CHECK_OUTPUT_CORRECT_TYPE_CALLED)
            self.assertTrue(not CHECK_CAN_CREATE_OUTPUT_CALLED)
            self.assertTrue(CHECK_OVERWRITE_EXISTING_FILES_CALLED)

def mock_check_input_exists(foo):
    global CHECK_INPUT_EXISTS_CALLED
    CHECK_INPUT_EXISTS_CALLED = True

def mock_check_input_readable(foo):
    global CHECK_INPUT_READABLE_CALLED
    CHECK_INPUT_READABLE_CALLED = True

def mock_check_input_correct_type(foo, bar):
    global CHECK_INPUT_CORRECT_TYPE_CALLED
    CHECK_INPUT_CORRECT_TYPE_CALLED = True

def mock_check_output_exists(foo, bar):
    global CHECK_OUTPUT_EXISTS_CALLED
    CHECK_OUTPUT_EXISTS_CALLED = True
    if not os.path.exists(foo):
        mock_check_can_create_output(foo)
    else:
        mock_check_output_correct_type(foo, bar)

def mock_check_can_create_output(foo):
    global CHECK_CAN_CREATE_OUTPUT_CALLED
    CHECK_CAN_CREATE_OUTPUT_CALLED = True

def mock_check_output_correct_type(foo, bar):
    global CHECK_OUTPUT_CORRECT_TYPE_CALLED
    CHECK_OUTPUT_CORRECT_TYPE_CALLED = True

def mock_check_tmpdir_exists(foo):
    global CHECK_TMPDIR_EXISTS_CALLED
    CHECK_TMPDIR_EXISTS_CALLED = True

def mock_check_overwrite_existing_files(foo, bar, baz, qux):
    global CHECK_OVERWRITE_EXISTING_FILES_CALLED
    CHECK_OVERWRITE_EXISTING_FILES_CALLED = True

def reset_called_methods():
    global CHECK_INPUT_EXISTS_CALLED
    global CHECK_INPUT_READABLE_CALLED
    global CHECK_INPUT_CORRECT_TYPE_CALLED
    global CHECK_OUTPUT_EXISTS_CALLED
    global CHECK_CAN_CREATE_OUTPUT_CALLED
    global CHECK_OUTPUT_CORRECT_TYPE_CALLED
    global CHECK_TMPDIR_EXISTS_CALLED
    global CHECK_OVERWRITE_EXISTING_FILES_CALLED

    CHECK_INPUT_EXISTS_CALLED = False
    CHECK_INPUT_READABLE_CALLED = False
    CHECK_INPUT_CORRECT_TYPE_CALLED = False
    CHECK_OUTPUT_EXISTS_CALLED = False
    CHECK_CAN_CREATE_OUTPUT_CALLED = False
    CHECK_OUTPUT_CORRECT_TYPE_CALLED = False
    CHECK_TMPDIR_EXISTS_CALLED = False
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

