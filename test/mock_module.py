#pylint: disable=invalid-name, global-statement, unused-argument
execute_called = False
report_called = False
validate_args_called = False
my_exception_string = None
predicted_output = "foo"

def init_mock():
    global execute_called
    execute_called = False
    global report_called
    report_called = False
    global validate_args_called
    validate_args_called = False
    global my_exception_string
    my_exception_string = None
    global predicted_output
    predicted_output = "foo"


def add_subparser(subparser):
    parser = subparser.add_parser("mock_module", help="foo")
    parser.add_argument("input", help="foo")
    parser.add_argument("output", help="foo")
    parser.add_argument("--force", action='store_true', help="foo")

def execute(args, execution_context):
    global execute_called
    execute_called = True
    if my_exception_string:
        raise Exception(my_exception_string)

def report_prediction(args):
    global report_called
    report_called = True
    if my_exception_string:
        raise Exception(my_exception_string)
    return predicted_output

def get_required_input_output_types():
    return ("directory", "directory")

def validate_args(args):
    global validate_args_called
    validate_args_called = True

init_mock()

