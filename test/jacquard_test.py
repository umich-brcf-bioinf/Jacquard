# pylint: disable=line-too-long, global-statement, unused-argument
# pylint: disable=invalid-name, too-many-locals, too-many-public-methods
# pylint: disable=too-few-public-methods
from __future__ import print_function, absolute_import, division

from argparse import Namespace
import os
import signal

from testfixtures import TempDirectory

import jacquard.jacquard as jacquard
import jacquard.utils.logger as logger
import jacquard.utils.utils as utils
import test.utils.mock_module as mock_module
import test.utils.test_case as test_case


MOCK_MOVE_TEMP_CONTENTS_CALLED = False
MOCK_PREFLIGHT_CALLED = False

def mock_preflight(output, desired_output_files, command):
    global MOCK_PREFLIGHT_CALLED
    MOCK_PREFLIGHT_CALLED = True

def mock_move_tmp_contents_to_original(args):
    global MOCK_MOVE_TEMP_CONTENTS_CALLED
    MOCK_MOVE_TEMP_CONTENTS_CALLED = True

def _change_mock_methods():
    jacquard._move_tmp_contents_to_original = mock_move_tmp_contents_to_original
    jacquard._preflight = mock_preflight


class MockSignalDispatcher(object):
    def __init__(self):
        self.calls = []

    def signal(self, signal_num, handler):
        self.calls.append((signal_num, handler))

class JacquardArgumentParserTestCase(test_case.JacquardBaseTestCase):
    def test_error_raisesMessage(self):
        message = "foo"
        parser = jacquard._JacquardArgumentParser()
        self.assertRaisesRegexp(utils.UsageError, "^foo$", parser.error, message)

    def test_error_raisesTransformedMessage(self):
        message = "argument subparser_name: invalid choice: 'foo' (choose from 'bar', 'baz')"
        parser = jacquard._JacquardArgumentParser()
        self.assertRaisesRegexp(utils.UsageError,
                                "^'foo' is not a Jacquard command; choose from 'bar', 'baz'.$",
                                parser.error,
                                message)

class JacquardTestCase(test_case.JacquardBaseTestCase):
    def test_version(self):
        self.assertEquals("1.1.2", jacquard.__version__)

    def test_get_execution_context(self):
        command = "foo input_dir output_dir"
        actual = jacquard._get_execution_context(command)
        expected = r'##jacquard=<timestamp=".*",command="foo input_dir output_dir",cwd=".*",source="Jacquard",version=".*">'
        self.assertRegexpMatches(actual[0], expected)

    def test_set_interrupt_handler(self):
        mock_signal_dispatcher = MockSignalDispatcher()

        jacquard._set_interrupt_handler(mock_signal_dispatcher.signal)

        self.assertEquals(2, len(mock_signal_dispatcher.calls))
        signal_type1, handler1 = mock_signal_dispatcher.calls[0]
        signal_type2, handler2 = mock_signal_dispatcher.calls[1]
        self.assertEquals(set([signal.SIGINT, signal.SIGTERM]),
                          set([signal_type1, signal_type2]))
        self.assertEquals(handler1, handler2)
        with self.assertRaises(SystemExit) as exit_code:
            handler1(None, None)
        self.assertEquals(exit_code.exception.code, 1)
        self.assertRegexpMatches(self.output.getvalue(), "WARNING: Jacquard was interrupted")

    def test_gracefulErrorMessageWhenUnanticipatedProblem(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = "I'm feeling angry"

            with self.assertRaises(SystemExit) as exit_code:
                jacquard._dispatch([mock_module], ["mock_module",
                                                   input_dir.path,
                                                   output_dir.path])
            self.assertEqual(1, exit_code.exception.code)

            actual_messages = self.output.getvalue().rstrip().split("\n")

            self.assertEquals(2, len(actual_messages))
            self.assertRegexpMatches(actual_messages[0], "I'm feeling angry")
            self.assertRegexpMatches(actual_messages[1], "Jacquard encountered an unanticipated problem.")

    def test_move_tmp_contents_to_original(self):
        with TempDirectory() as parent_dir:
            output_dir = os.path.join(parent_dir.path, "my_data")
            tmp_dir = os.path.join(parent_dir.path, "tmp")
            args = Namespace(original_output=output_dir,
                             temp_working_dir=tmp_dir)

            os.mkdir(tmp_dir)
            os.mkdir(os.path.join(tmp_dir, "my_data"))
            with open(os.path.join(tmp_dir, "my_data", "A.txt"), "w") as file_a, \
                    open(os.path.join(tmp_dir, "my_data", "B.txt"), "w") as file_b:
                file_a.write("A")
                file_b.write("B")

            jacquard._move_tmp_contents_to_original(args)

            actual_files = os.listdir(output_dir)
            self.assertEquals(2, len(actual_files))
            self.assertEquals(["A.txt", "B.txt"], sorted(actual_files))
            self.assertFalse(os.path.isdir(tmp_dir))

    def test_move_tmp_contents_to_original_mergesToExistingDir(self):
        with TempDirectory() as parent_dir:
            dest_dir = os.path.join(parent_dir.path, "my_data")
            tmp_dir = os.path.join(parent_dir.path, "tmp")
            args = Namespace(original_output=dest_dir,
                             temp_working_dir=tmp_dir)

            os.mkdir(tmp_dir)
            os.mkdir(dest_dir)

            os.mkdir(os.path.join(tmp_dir, "my_data"))
            with open(os.path.join(tmp_dir, "my_data", "A.txt"), "w") as file_a, \
                    open(os.path.join(tmp_dir, "my_data", "B.txt"), "w") as file_b, \
                    open(os.path.join(dest_dir, "C.txt"), "w") as file_c:
                file_a.write("A")
                file_b.write("B")
                file_c.write("C")

            jacquard._move_tmp_contents_to_original(args)

            actual_files = os.listdir(dest_dir)
            self.assertEquals(3, len(actual_files))
            self.assertEquals(["A.txt", "B.txt", "C.txt"], sorted(actual_files))
            self.assertFalse(os.path.isdir(tmp_dir))


class JacquardTestCase_dispatchOnly(test_case.JacquardBaseTestCase):
    def setUp(self):
        test_case.JacquardBaseTestCase.setUp(self)
        self.original_move_tmp_contents = jacquard._move_tmp_contents_to_original
        _change_mock_methods()

    def tearDown(self):
        jacquard._move_tmp_contents_to_original = self.original_move_tmp_contents
        test_case.JacquardBaseTestCase.tearDown(self)

    def test_dispatch(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            jacquard._dispatch([mock_module], ["mock_module",
                                               input_dir.path,
                                               output_dir.path])
            self.assertTrue(mock_module.execute_called)
            self.assertTrue(mock_module.report_called)
            self.assertTrue(MOCK_MOVE_TEMP_CONTENTS_CALLED)

    def test_dispatch_unknownCommand(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            with self.assertRaises(SystemExit) as exit_code:
                jacquard._dispatch([mock_module], ["foo_command",
                                                   input_dir.path,
                                                   output_dir.path])
            self.assertEquals(1, exit_code.exception.code)
            self.assertRegexpMatches(self.output.getvalue(),
                                     r"'foo_command' is not a Jacquard command; choose from 'mock_module'")


    def test_dispatch_nonEmptyOutputDir(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            output_dir.write("file1.vcf", b"foo")
            output_file = os.path.join(output_dir.path, "file1.vcf")
            mock_module.my_exception_string = ""
            mock_module.predicted_output = set(["file1.vcf"])
            with self.assertRaises(SystemExit) as exit_code:
                jacquard._dispatch([mock_module], ["mock_module",
                                                   input_dir.path,
                                                   output_file])

            self.assertEqual(1, exit_code.exception.code)

    def test_dispatch_forceNonEmptyOutputDir(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            output_dir.write("file1.vcf", b"foo")
            mock_module.my_exception_string = ""

            jacquard._dispatch([mock_module], ["mock_module",
                                               input_dir.path,
                                               output_dir.path,
                                               "--force"])

            self.assertTrue(1 == 1, "Force does not result in premature exit.")

    def test_dispatch_done(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            jacquard._dispatch([mock_module], ["mock_module",
                                               input_dir.path,
                                               output_dir.path])
            actual_messages = self.output.getvalue().rstrip().split("\n")
            self.assertRegexpMatches(actual_messages[-1], "Done")

    def test_dispatch_doneWithWarnings(self):
        with TempDirectory() as input_dir, TempDirectory() as output_dir:
            mock_module.my_exception_string = ""
            logger.WARNING_OCCURRED = True
            jacquard._dispatch([mock_module], ["mock_module",
                                               input_dir.path,
                                               output_dir.path])
            actual_messages = self.output.getvalue().rstrip().split("\n")
            self.assertRegexpMatches(actual_messages[-1], r"Done. \(See warnings above\)")
            logger.WARNING_OCCURRED = False

class JacquardFunctionalTestCase(test_case.JacquardBaseTestCase):
    def test_functional_jacquard(self):
        with TempDirectory() as output_dir:
            file_dirname = os.path.dirname(os.path.realpath(__file__))
            module_testdir = os.path.join(file_dirname,
                                          "functional_tests",
                                          "jacquard_test")

            initial_input = os.path.join(module_testdir, "00_input")
            outputs = {"01_translate": "translated",
                       "02_merge": "merged.vcf",
                       "03_summarize": "summarized.vcf",
                       "04_expand": "expanded.tsv"}

            translate_out = os.path.join(output_dir.path, outputs["01_translate"])
            merge_out = os.path.join(output_dir.path, outputs["02_merge"])
            summarize_out = os.path.join(output_dir.path, outputs["03_summarize"])
            expanded_out = os.path.join(output_dir.path, outputs["04_expand"])

            commands = [["translate", initial_input, translate_out, "--force"],
                        ["merge", translate_out, merge_out, "--force"],
                        ["summarize", merge_out, summarize_out, "--force"],
                        ["expand", summarize_out, expanded_out, "--force"]]

            prefixes = ["01_", "02_", "03_", "04_"]
            count = 0
            for command in commands:
                command_name = prefixes[count] + command[0]
                expected_dir = os.path.join(module_testdir,
                                            command_name,
                                            "benchmark",
                                            outputs[command_name])
                self.assertCommand(command, expected_dir)

                count += 1
