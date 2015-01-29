from __future__ import absolute_import

import glob
import os
import errno

import jacquard.utils as utils

TMP_DIR_NAME = "jacquard_tmp"


def preflight(args, module):
    input_path = args.input
    output_path = args.output
    module_name = args.subparser_name
    force_flag = args.force

    required_input_type, required_output_type =\
                                    module.get_required_input_output_types()

    _check_input_exists(input_path)
    _check_input_readable(input_path)
    _check_input_correct_type(input_path, required_input_type)
    _check_output_exists(output_path, required_output_type)
    _check_tmpdir_exists(output_path)

    predicted_output = module.report_prediction(args)
    _check_overwrite_existing_files(output_path,
                                    predicted_output,
                                    module_name,
                                    force_flag)

def _makepath(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def _check_input_exists(input_path):
    if not os.path.exists(input_path):
        raise utils.JQException(("Specified input [{}] does not exist. Review "\
                                 "inputs and try again.").format(input_path))

def _check_input_readable(input_path):
    try:
        if os.path.isdir(input_path):
            os.listdir(input_path)
        else:
            open(input_path, "r").close()
    except (OSError, IOError):
        raise utils.JQException(("Specified input [{}] cannot be read. Review "\
                                 "inputs and try again.").format(input_path))

def _check_input_correct_type(input_path, required_type):
    if (required_type == "file" and not os.path.isfile(input_path)) or\
    (required_type == "directory" and not os.path.isdir(input_path)):
        raise utils.JQException(("Specified input [{}] does not match "\
                                 "command's required file type. "\
                                 "Review inputs and try again.")\
                                 .format(input_path))

def _check_output_exists(output_path, required_type):
    if not os.path.exists(output_path):
        _check_can_create_output(output_path)
    else:
        _check_output_correct_type(output_path, required_type)

def _check_can_create_output(output_path):
    output_path = os.path.dirname(output_path)
    try:
        _makepath(output_path)
    except OSError:
        raise utils.JQException(("Specified output [{}] does not exist "
                                 "and cannot be created. Review inputs "
                                 "and try again.").format(output_path))

def _check_output_correct_type(output_path, required_type):
    if (required_type == "file" and os.path.isdir(output_path)) or\
    (required_type == "directory" and not os.path.isdir(output_path)):
        raise utils.JQException(("Specified output [{}] does not match "
                                 "command's required file type. "
                                 "Review inputs and try again.")\
                                 .format(output_path))

def _check_tmpdir_exists(output_path):
    if os.path.isdir(output_path):
        output_dir = output_path
    else:
        output_dir = os.path.dirname(output_path)

    try:
        os.listdir(output_dir)
    except OSError:
        raise utils.JQException(("Output directory [{}] cannot be read."\
                                 "Review inputs and try again.")\
                                 .format(output_dir))

    if TMP_DIR_NAME not in os.listdir(output_dir):
        try:
            tmp_output = os.path.join(output_dir, TMP_DIR_NAME)
            os.mkdir(tmp_output)
        except OSError:
            raise utils.JQException(("A tmp directory does not exist in "\
                                     "[{}] and cannot be created. Review "\
                                     "inputs and try again.")\
                                     .format(output_path))

def _check_overwrite_existing_files(output, predicted_output, command, force=0):
    if os.path.isdir(output):
        existing_output_paths = sorted(glob.glob(os.path.join(output, "*.vcf")))
    else:
        existing_output_paths = [output]
    existing_output = set([os.path.basename(i) for i in existing_output_paths])

    intersection = existing_output.intersection(predicted_output)
    if intersection and not force:
        raise utils.JQException(("ERROR: The command [{}] would "
                                "overwrite existing files {}; review "
                                "command/output dir to avoid overwriting or "
                                "use the flag '--force'. Type 'jacquard -h' "
                                "for more details").format(command,
                                                           list(intersection)))
