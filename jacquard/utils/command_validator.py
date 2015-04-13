"""Validates command preconditions.

Specifically checks that the command, arguments, and environment
(e.g. input/output directories or files) are consistent and plausible.

Each validation function evaluates a specific precondition.
Each function is allowed to:
 * change the environment (e.g. create a dir)
 * change the args (replace the original output dir with a new temp output dir)
 * add arguments which may be required for sub-commands
 * delegate to/interact with a sub-command
 * raise a UsageException if things look problematic
"""
from __future__ import print_function, absolute_import, division

import errno
import glob
import os
import time

import jacquard.utils.utils as utils


_TEMP_WORKING_DIR_FORMAT = "jacquard.{}.{}.tmp"


def _actual_type(path):
    if os.path.isdir(path):
        return "directory"
    else:
        return "file"

def _build_collision_message(command, collisions):
    total_collisions = len(collisions)
    if total_collisions == 1:
        return ("The {} command would "
                "overwrite the existing file [{}]; review "
                "command/output dir to avoid overwriting or "
                "use the flag '--force'.").format(command,
                                                  collisions[0])
    cutoff = 5
    collision_list = ", ".join(collisions[0:min(cutoff, total_collisions)])
    if total_collisions > cutoff:
        omitted = total_collisions - cutoff
        collision_list += ", ...({} file(s) omitted)".format(omitted)
    return ("The {} command would "
            "overwrite {} existing files [{}]; review "
            "command/output dir to avoid overwriting or "
            "use the flag '--force'.").format(command,
                                              total_collisions,
                                              collision_list)

def _check_input_correct_type(dummy, args):
    module_name = args.subparser_name
    input_path = args.input
    required_type = args.required_input_type
    actual_type = _actual_type(input_path)
    if required_type != actual_type:
        raise utils.UsageError(("The {} command requires a {} as "
                                "input, but the specified input [{}] is a {}. "
                                "Review inputs and try again.") \
                                 .format(module_name,
                                         required_type,
                                         input_path,
                                         actual_type))

def _check_input_exists(dummy, args):
    if not os.path.exists(args.input):
        raise utils.UsageError(("Specified input [{}] does not exist. Review "\
                                 "inputs and try again.").format(args.input))

def _check_input_readable(dummy, args):
    try:
        if os.path.isdir(args.input):
            os.listdir(args.input)
        else:
            open(args.input, "r").close()
    except (OSError, IOError):
        raise utils.UsageError(("Specified input [{}] cannot be read. Review "
                                "inputs and try again.").format(args.input))

def _check_output_correct_type(module_name, output_path, required_type):
    actual_type = _actual_type(output_path)
    if required_type != actual_type:
        raise utils.UsageError(("The {} command outputs a {}, but the "
                                "specified output [{}] is a {}. "
                                "Review inputs and try again.")\
                                 .format(module_name,
                                         required_type,
                                         output_path,
                                         actual_type))

def _check_output_exists(dummy, args):
    if os.path.exists(args.output_path):
        _check_output_correct_type(args.subparser_name,
                                   args.output_path,
                                   args.required_output_type)

def _check_overwrite_existing_files(module, args):
    output = args.output_path
    if not os.path.isdir(output):
        output = os.path.dirname(output)

    existing_output_paths = sorted(glob.glob(os.path.join(output, "*")))
    existing_output = set([os.path.basename(i) for i in existing_output_paths])
    predicted_output = module.report_prediction(args)
    collisions = sorted(list(existing_output.intersection(predicted_output)))

    if collisions and not args.force:
        message = _build_collision_message(args.subparser_name, collisions)
        raise utils.UsageError(message)

def _check_there_will_be_output(module, args):
    predicted_output = module.report_prediction(args)

    if not predicted_output:
        message = ("Executing the {} command with the input [{}] would not "
                   "create any output files. Review inputs and try again.")\
                   .format(args.subparser_name, args.input)
        raise utils.UsageError(message)

def _check_valid_args(module, args):
    module.validate_args(args)

def _create_temp_working_dir(dummy, args):
    try:
        _makepath(args.temp_working_dir)
        if args.required_output_type == "directory":
            _makepath(args.output)
    except OSError:
        parent_dir = os.path.dirname(args.temp_working_dir)
        raise utils.UsageError(("Jacquard cannot write to output directory "
                                "[{}]. Review inputs and try again.")\
                                .format(parent_dir))

def _set_temp_working_dir(dummy, args):
    original_output = args.original_output
    required_output = args.required_output_type
    abs_original_output = os.path.abspath(original_output)
    pid = os.getpid()
    microseconds_since_epoch = int(time.time() * 1000 * 1000)
    dir_name = _TEMP_WORKING_DIR_FORMAT.format(str(pid),
                                               str(microseconds_since_epoch))

    place_temp_inside_output = True
    if required_output == "file":
        place_temp_inside_output = False
    elif required_output == "directory" and not os.path.isdir(original_output):
        place_temp_inside_output = False

    if place_temp_inside_output:
        base_dir = abs_original_output
        temp_working_dir = os.path.join(base_dir, dir_name)
        new_output = temp_working_dir
    else:
        base_dir = os.path.dirname(abs_original_output)
        temp_working_dir = os.path.join(base_dir, dir_name)
        new_output = os.path.join(temp_working_dir,
                                  os.path.basename(abs_original_output))

    args.temp_working_dir = temp_working_dir
    args.output = new_output

def _makepath(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def _set_required_types(module, args):
    (args.required_input_type,
     args.required_output_type) = module.get_required_input_output_types()

def _set_output_paths(dummy, args):
    args.original_output = args.output
    args.output_path = os.path.abspath(args.original_output)

_VALIDATION_TASKS = [_set_output_paths,
                     _set_required_types,
                     _set_temp_working_dir,
                     _check_input_exists,
                     _check_input_readable,
                     _check_input_correct_type,
                     _check_output_exists,
                     _create_temp_working_dir,
                     _check_there_will_be_output,
                     _check_overwrite_existing_files,
                     _check_valid_args]

def preflight(command, args):
    for validate in _VALIDATION_TASKS:
        validate(command, args)
