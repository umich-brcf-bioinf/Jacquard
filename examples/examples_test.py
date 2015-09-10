#pylint: disable=too-many-public-methods, invalid-name
import os
import zipfile

from testfixtures import TempDirectory

import test.utils.test_case as test_case

def unzip_directory(output_dir):
    file_dirname = os.path.dirname(os.path.realpath(__file__))
    zipped_dirname = os.path.join(file_dirname, "examples.zip")
    with zipfile.ZipFile(zipped_dirname, "r") as z:
        z.extractall(output_dir)


class ExamplesFunctionalTest(test_case.JacquardBaseTestCase):
    def test_examples(self):
        with TempDirectory() as output_dir, TempDirectory() as reference_dir:
            unzip_directory(reference_dir.path)

            initial_input = os.path.join(reference_dir.path,
                                         "examples",
                                        "00-input_vcfs")
            outputs = {"translate": "01-translated",
                       "merge": "02-merged.vcf",
                       "summarize": "03-summarized.vcf",
                       "expand": "04-expanded.tsv"}

            translate_out = os.path.join(output_dir.path, outputs["translate"])
            merge_out = os.path.join(output_dir.path, outputs["merge"])
            summarize_out = os.path.join(output_dir.path, outputs["summarize"])
            expanded_out = os.path.join(output_dir.path, outputs["expand"])

            commands = [["translate", initial_input, translate_out, "--force"],
                        ["merge", translate_out, merge_out, "--force"],
                        ["summarize", merge_out, summarize_out, "--force"],
                        ["expand", summarize_out, expanded_out, "--force"]]

            for command in commands:
                command_name = command[0]
                expected_dir = os.path.join(reference_dir.path,
                                            "examples",
                                            outputs[command_name])
                self.assertCommand(command, expected_dir)

