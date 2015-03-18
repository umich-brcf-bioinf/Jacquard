#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=global-statement, unused-argument, too-few-public-methods
#pylint: disable=too-many-instance-attributes, no-member
from __future__ import absolute_import

from argparse import Namespace
import os
import re
import shutil
import subprocess
import sys
import time

from testfixtures import TempDirectory

import jacquard.command_validator as command_validator
import jacquard.utils as utils
import test.mock_module as mock_module
import test.test_case as test_case


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

class CommandValidatorTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(CommandValidatorTestCase, self).setUp()
        mock_module.init_mock()

    def tearDown(self):
        super(CommandValidatorTestCase, self).tearDown()

    def test_validation_tasks(self):
        self.assertEquals(12, len(command_validator._VALIDATION_TASKS))

    def test_check_output_exists_fileExistsCorrectType(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            args = Namespace(subparser_name="awesomeCommand",
                             output_path=output_file,
                             required_output_type="file")
            command_validator._check_output_exists(None, args)
            self.ok()

    def test_check_output_exists_fileExistsIncorrectType(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file.vcf", "foo")
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            args = Namespace(subparser_name="awesomeCommand",
                             output_path=output_file,
                             required_output_type="directory")
            self.assertRaisesRegexp(utils.UsageError,
                                    "The awesomeCommand command outputs a directory, but the specified output.*output_file.vcf.*is a file",
                                    command_validator._check_output_exists,
                                    None,
                                    args)

    def test_check_output_exists_fileDoesNotExist(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path, "output_file.vcf")
            args = Namespace(subparser_name="awesomeCommand",
                             output_path=output_file,
                             required_output_type="x")
            command_validator._check_output_exists(None, args)
            self.ok()

    def test_check_output_exists_directoryExistsCorrectType(self):
        with TempDirectory() as output_dir:
            output_dir.makedir("output_dir")
            expected_output_dir = os.path.join(output_dir.path, "output_dir")
            args = Namespace(subparser_name="awesomeCommand",
                             output_path=expected_output_dir,
                             required_output_type="directory")
            command_validator._check_output_exists(None, args)
            self.ok()

    def test_check_output_exists_directoryExistsIncorrectType(self):
        with TempDirectory() as output_dir:
            args = Namespace(subparser_name="awesomeCommand",
                             output_path=output_dir.path,
                             required_output_type="file")
            self.assertRaisesRegexp(utils.UsageError,
                                    "The awesomeCommand command outputs a file, but the specified output.*is a directory",
                                    command_validator._check_output_exists,
                                    None,
                                    args)

    def test_set_output_paths(self):
        with TempDirectory() as output_dir:
            args = Namespace(output=output_dir.path)
            command_validator._set_output_paths(None, args)
            self.assertEquals(args.original_output, output_dir.path)
            self.assertEquals(args.output_path, os.path.abspath(output_dir.path))

    def test_get_temp_working_dir_insideDestIfDestDirExists(self):
        with TempDirectory() as original_output_dir:
            original_output = os.path.join(original_output_dir.path, "dest_dir")
            os.mkdir(original_output)
            args = Namespace(original_output=original_output,
                             required_output_type="directory")

            command_validator._set_temp_working_dir(None, args)
        self.assertIn(original_output, args.temp_working_dir)
        self.assertEquals(args.temp_working_dir, args.output)

    def test_set_temp_working_dir_besideDestIfDestDirNotExists(self):
        with TempDirectory() as original_parent_dir:
            parent_dir = original_parent_dir.path
            args = Namespace(original_output=os.path.join(parent_dir, "dest_dir"),
                             required_output_type="directory")

            command_validator._set_temp_working_dir(None, args)

        self.assertIn(parent_dir, args.temp_working_dir)
        self.assertIn(args.temp_working_dir, args.output)

    def test_get_temp_working_dir_besideDestIfDestIsFile(self):
        with TempDirectory() as original_parent_dir:
            parent_dir = original_parent_dir.path
            original_output = os.path.join(parent_dir, "dest_file")
            args = Namespace(original_output=original_output,
                             required_output_type="file")
            command_validator._set_temp_working_dir(None, args)

        self.assertIn(parent_dir, args.temp_working_dir)
        self.assertIn(args.temp_working_dir, args.output)

    def test_get_temp_working_dir_distinct(self):
        orig_output_path = os.path.join("foo", "bar")
        args1 = Namespace(original_output=orig_output_path,
                          required_output_type="directory")
        command_validator._set_temp_working_dir(None, args1)
        time.sleep(0.01)
        args2 = Namespace(original_output=orig_output_path,
                          required_output_type="directory")
        command_validator._set_temp_working_dir(None, args2)
        self.assertNotEqual(args1.temp_working_dir, args2.temp_working_dir)

    def test_get_temp_working_dir_complexPathOk(self):
        original_output = os.path.join("foo.bar", "baz.hoopy", "frood")
        args = Namespace(original_output=original_output,
                         required_output_type="directory")
        command_validator._set_temp_working_dir(None, args)
        self.assertIn(os.path.dirname(original_output), args.temp_working_dir)
        self.assertIn(args.temp_working_dir, args.output)

    def test_get_temp_working_dir_trailingSlashOk(self):
        args = Namespace(original_output=os.path.join("foo", "bar") + os.sep,
                         required_output_type="directory")
        command_validator._set_temp_working_dir(None, args)
        self.assertRegexpMatches(args.temp_working_dir,
                                 r"foo.*jacquard.*tmp")
        self.assertRegexpMatches(args.output,
                                 "foo.*jacquard.*tmp.*bar")

    def test_get_temp_working_dir_currentDirOk(self):
        original_cwd = os.getcwd()
        try:
            with TempDirectory() as cwd:
                cwd_absolute_path = os.path.abspath(cwd.path)
                os.chdir(cwd_absolute_path)
                args = Namespace(original_output=os.path.join("."),
                                 required_output_type="directory")

                command_validator._set_temp_working_dir(None, args)
                os.chdir(original_cwd)
                self.assertIn(cwd_absolute_path, args.temp_working_dir)
                self.assertEquals(args.temp_working_dir, args.output)
        finally:
            os.chdir(original_cwd)

    def test_get_temp_working_dir_relativePathMadeAbsolute(self):
        original_cwd = os.getcwd()
        try:
            with TempDirectory() as cwd:
                cwd_absolute_path = os.path.abspath(cwd.path)
                os.chdir(cwd_absolute_path)
                args = Namespace(original_output=os.path.join(".", "foo", "bar"),
                                 required_output_type="directory")
                command_validator._set_temp_working_dir(None, args)
                os.chdir(original_cwd)
                self.assertIn(cwd_absolute_path, args.temp_working_dir)
        finally:
            os.chdir(original_cwd)

    def test_check_input_exists_fileExists(self):
        with TempDirectory() as input_dir:
            input_dir.write("existing_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "existing_file.vcf")
            command_validator._check_input_exists(None,
                                                  Namespace(input=input_file))
            self.assertTrue(1 == 1)

    def test_check_input_exists_fileNonExists(self):
        with TempDirectory() as input_dir:
            input_file = os.path.join(input_dir.path, "nonexistant_file.vcf")
            self.assertRaisesRegexp(utils.UsageError,
                                    "Specified input .* does not exist",
                                    command_validator._check_input_exists,
                                    None,
                                    Namespace(input=input_file))

    def test_check_input_exists_directoryExists(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("existing_dir")
            input_file = os.path.join(input_dir.path, "existing_dir")
            command_validator._check_input_exists(None,
                                                  Namespace(input=input_file))
            self.assertTrue(1 == 1)

    def test_check_input_exists_directoryNonExists(self):
        with TempDirectory() as input_dir:
            input_file = os.path.join(input_dir.path, "missing_dir")
            self.assertRaisesRegexp(utils.UsageError,
                                    "Specified input .* does not exist",
                                    command_validator._check_input_exists,
                                    None,
                                    Namespace(input=input_file))

    def test_check_input_readable_fileReadable(self):
        with TempDirectory() as input_dir:
            input_dir.write("readable_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "readable_file.vcf")
            command_validator._check_input_readable(None,
                                                    Namespace(input=input_file))

            self.assertTrue(1 == 1)

    def test_check_input_readable_fileNotReadable(self):
        with TempDirectory() as input_dir:
            input_dir.write("readable_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "readable_file.vcf")

            try:
                make_unreadable(input_file)
                self.assertRaisesRegexp(utils.UsageError,
                                        "Specified input .* cannot be read",
                                        command_validator._check_input_readable,
                                        None,
                                        Namespace(input=input_file))

            finally:
                cleanup_unreadable(input_file)

    def test_check_input_readable_directoryReadable(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("readable_directory")
            readable_dir = os.path.join(input_dir.path, "readable_directory")
            command_validator._check_input_readable(None,
                                                    Namespace(input=readable_dir))

            self.assertTrue(1 == 1)

    def test_check_input_readable_directoryNotReadable(self):
        with TempDirectory() as input_dir:
            input_dir.makedir("unreadable_directory")
            unreadable_dir = os.path.join(input_dir.path, "unreadable_directory")

            try:
                make_unreadable(unreadable_dir)
                self.assertRaisesRegexp(utils.UsageError,
                                        "Specified input .* cannot be read",
                                        command_validator._check_input_readable,
                                        None,
                                        Namespace(input=unreadable_dir))

            finally:
                cleanup_unreadable(unreadable_dir)

    def test_check_input_correct_type_correct(self):
        with TempDirectory() as input_dir:
            args_requires_dir = Namespace(subparser_name="awesomeCommand",
                                          input=input_dir.path,
                                          required_input_type="directory")
            command_validator._check_input_correct_type(None,
                                                        args_requires_dir)
            self.assertTrue(1 == 1)

            input_dir.write("input_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "input_file.vcf")
            args_requires_file = Namespace(subparser_name="awesomeCommand",
                                           input=input_file,
                                           required_input_type="file")
            command_validator._check_input_correct_type(None,
                                                        args_requires_file)
            self.assertTrue(1 == 1)

    def test_check_input_correct_type_incorrect(self):
        with TempDirectory() as input_dir:
            input_dir.write("input_file.vcf", "foo")
            input_file = os.path.join(input_dir.path, "input_file.vcf")
            args_requires_dir = Namespace(subparser_name="awesomeCommand",
                                          input=input_file,
                                          required_input_type="directory")
            self.assertRaisesRegexp(utils.UsageError,
                                    r"The awesomeCommand command requires a directory.*input_file.vcf.*is a file",
                                    command_validator._check_input_correct_type,
                                    None,
                                    args_requires_dir)

            args_required_file = Namespace(subparser_name="awesomeCommand",
                                           input=input_dir.path,
                                           required_input_type="file")
            self.assertRaisesRegexp(utils.UsageError,
                                    r"The awesomeCommand command requires a file.*"+re.escape(input_dir.path)+r".*is a directory.*",
                                    command_validator._check_input_correct_type,
                                    None,
                                    args_required_file)

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
                                    "The awesomeCommand command outputs a file.*" + re.escape(output_dir.path) + ".*is a directory",
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

            command_validator._create_temp_working_dir(None, args)
            self.assertTrue(1 == 1)

    def test_create_temp_working_dir_createsIfNotExists(self):
        with TempDirectory() as output_dir:
            tmp_dir = os.path.join(output_dir.path, "jacquard.tmp")
            args = Namespace(output=os.path.join(output_dir.path, "dest_dir"),
                             temp_working_dir=tmp_dir,
                             required_output_type="directory")

            command_validator._create_temp_working_dir(None, args)
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

            command_validator._create_temp_working_dir(None, args)
            self.ok()

    def test_check_overwrite_existing_files_whenNotOverwriteWillNotRaise(self):
        with TempDirectory() as output_dir:
            args = Namespace(output_path=output_dir.path,
                             subparser_name="awesomeCommand")
            mock_module.predicted_output = "foo.vcf"
            command_validator._check_overwrite_existing_files(mock_module,
                                                              args)
            self.ok()

    def test_check_overwrite_existing_files_whenOverwriteSingleFileWillRaise(self):
        with TempDirectory() as output_dir:
            output_file = os.path.join(output_dir.path, "output_file")
            output_dir.write("output_file", "foo")
            args = Namespace(output_path=output_file,
                             subparser_name="awesomeCommand",
                             force=0)
            mock_module.predicted_output = set([os.path.basename(output_file)])
            self.assertRaisesRegexp(utils.UsageError,
                                    ("The awesomeCommand command would "
                                     "overwrite the existing file "
                                     r"\[output_file\]; review "
                                     "command/output dir to avoid overwriting or "
                                     "use the flag '--force'."),
                                    command_validator._check_overwrite_existing_files,
                                    mock_module,
                                    args)

    def test_check_overwrite_whenOverwriteManyFilesWillRaise(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file1.vcf", "foo")
            output_dir.write("output_file2.vcf", "foo")
            args = Namespace(output_path=output_dir.path,
                             subparser_name="awesomeCommand",
                             force=0)
            mock_module.predicted_output = set(["output_file1.vcf",
                                                "output_file2.vcf"])
            self.assertRaisesRegexp(utils.UsageError,
                                    "The awesomeCommand command would overwrite 2 "
                                    r"existing files \[output_file1.vcf, output_file2.vcf\]; "
                                    "review command/output dir to avoid overwriting or "
                                    "use the flag '--force'.",
                                    command_validator._check_overwrite_existing_files,
                                    mock_module,
                                    args)

    def test_check_overwrite_existing_files_whenOverwriteMoreThanFiveFilesWillRaiseAbbreviatedUsageError(self):
        with TempDirectory() as output_dir:
            args = Namespace(output_path=output_dir.path,
                             subparser_name="awesomeCommand",
                             force=0)
            predicted_outputs = set()
            for x in xrange(1, 7):
                file_name = "out{}.vcf".format(str(x))
                output_dir.write(file_name, "foo")
                predicted_outputs.add(file_name)
            mock_module.predicted_output = predicted_outputs
            self.assertRaisesRegexp(utils.UsageError,
                                    "The awesomeCommand command would overwrite 6 "
                                    r"existing files \[out1.vcf, out2.vcf, out3.vcf, out4.vcf, out5.vcf, "
                                    r"...\(1 file\(s\) omitted\)\]; "
                                    "review command/output dir to avoid "
                                    "overwriting or use the flag '--force'.",
                                    command_validator._check_overwrite_existing_files,
                                    mock_module,
                                    args)

    def test_check_overwrite_existing_files_whenNotOverwriteDirectoryWillNotRaise(self):
        with TempDirectory() as output_dir:
            args = Namespace(output_path=output_dir.path,
                             subparser_name="awesomeCommand",
                             force=0)
            mock_module.predicted_output = set(["output_file1.vcf",
                                                "output_file2.vcf"])
            command_validator._check_overwrite_existing_files(mock_module,
                                                              args)
            self.ok()

    def test_check_overwrite_whenForcedOnExistingFileDoesNotRaise(self):
        with TempDirectory() as output_dir:
            output_dir.write("output_file1.vcf", "foo")
            output_dir.write("output_file2.vcf", "foo")
            args = Namespace(output_path=output_dir.path,
                             subparser_name="awesomeCommand",
                             force=1)
            mock_module.predicted_output = set(["output_file1.vcf",
                                                "output_file2.vcf"])
            command_validator._check_overwrite_existing_files(mock_module,
                                                              args)
            self.ok()

    def test_check_there_will_be_output(self):
        mock_module.predicted_output = ["out1", "out2"]
        args = Namespace(subparser_name="awesomeCommand", input="input_dir")
        command_validator._check_there_will_be_output(mock_module, args)
        self.ok()

    def test_check_there_will_be_output_raisesIfNoPredictedOutput(self):
        mock_module.predicted_output = []
        args = Namespace(subparser_name="awesomeCommand", input="input_dir")
        expected_regexp = ("Executing the awesomeCommand command with the "
                           r"input \[input_dir\] would not create any output "
                           "files. Review inputs and try again.")
        self.assertRaisesRegexp(utils.UsageError,
                                expected_regexp,
                                command_validator._check_there_will_be_output,
                                mock_module,
                                args)

    def test_check_valid_args(self):
        args = Namespace()
        command_validator._check_valid_args(mock_module, args)
        self.assertTrue(mock_module.validate_args_called)

    def test_check_snp_indel_pairing_missingIndel(self):
        with TempDirectory() as input_dir:
            input_dir.write("patientA.snp.vcf", "foo")
            input_dir.write("patientA.indel.vcf", "foo")
            input_dir.write("patientB.snp.vcf", "foo")
            args = Namespace(input=input_dir.path,
                             subparser_name="awesomeCommand",
                             force=0)

            self.assertRaisesRegexp(utils.UsageError,
                                    "Some VCFs were missing either a snp/snvs or an indel/indels file. Review inputs/command options and try again.",
                                    command_validator._check_input_snp_indel_pairing,
                                    None,
                                    args)

    def test_check_snp_indel_pairing_missingSnvs(self):
        with TempDirectory() as input_dir:
            input_dir.write("patientA.snp.vcf", "foo")
            input_dir.write("patientA.indel.vcf", "foo")
            input_dir.write("patientB.indels.vcf", "foo")
            args = Namespace(input=input_dir.path,
                             subparser_name="awesomeCommand",
                             force=0)

            self.assertRaisesRegexp(utils.UsageError,
                                    "Some VCFs were missing either a snp/snvs or an indel/indels file. Review inputs/command options and try again.",
                                    command_validator._check_input_snp_indel_pairing,
                                    None,
                                    args)

    def test_check_snp_indel_pairing_missingSnvsSamePatient(self):
        with TempDirectory() as input_dir:
            input_dir.write("patientA.snp.vcf", "foo")
            input_dir.write("patientA.indel.vcf", "foo")
            input_dir.write("patientA.indels.vcf", "foo")
            args = Namespace(input=input_dir.path,
                             subparser_name="awesomeCommand",
                             force=0)

            self.assertRaisesRegexp(utils.UsageError,
                                    "Some VCFs were missing either a snp/snvs or an indel/indels file. Review inputs/command options and try again.",
                                    command_validator._check_input_snp_indel_pairing,
                                    None,
                                    args)

    def test_check_snp_indel_pairing_allSnpOkay(self):
        with TempDirectory() as input_dir:
            input_dir.write("patientA.snp.vcf", "foo")
            input_dir.write("patientB.snp.vcf", "foo")
            args = Namespace(input=input_dir.path,
                             subparser_name="awesomeCommand",
                             force=0)

            command_validator._check_input_snp_indel_pairing(None, args)
            self.assertTrue(1, 1)


class MockTask(object):
    def __init__(self, error_message=None):
        self.error_message = error_message
        self.module = None
        self.args = None
        self.execute_called = False

    def execute(self, module, args):
        self.execute_called = True
        if self.error_message:
            raise Exception(self.error_message)
        else:
            self.module = module
            self.args = args

class CommandValidatorPreflightTestCase(test_case.JacquardBaseTestCase):
    def setUp(self):
        super(CommandValidatorPreflightTestCase, self).setUp()
        self.mock_module_original_predicted_output = mock_module.predicted_output
        self.command_validator_original_validation_tasks = list(command_validator._VALIDATION_TASKS)


    def tearDown(self):
        mock_module.predicted_output = self.mock_module_original_predicted_output
        command_validator._VALIDATION_TASKS = self.command_validator_original_validation_tasks
        super(CommandValidatorPreflightTestCase, self).tearDown()

    def test_preflight_runsAllValidationTasks(self):
        args = Namespace()
        task1 = MockTask()
        task2 = MockTask()
        command_validator._VALIDATION_TASKS = [task1.execute,
                                               task2.execute]
        command_validator.preflight(mock_module, args)

        self.assertTrue(task1.execute_called)
        self.assertEqual(task1.module, mock_module)
        self.assertEqual(task1.args, args)
        self.assertTrue(task2.execute_called)
        self.assertEqual(task2.module, mock_module)
        self.assertEqual(task2.args, args)

    def test_preflight_failsOnFirstError(self):
        args = Namespace()
        task1 = MockTask("error")
        task2 = MockTask()
        command_validator._VALIDATION_TASKS = [task1.execute,
                                               task2.execute]
        self.assertRaises(Exception,
                          command_validator.preflight,
                          mock_module,
                          args)

        self.assertTrue(task1.execute_called)
        self.assertFalse(task2.execute_called)


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

