#pylint: disable=invalid-name, global-statement, unused-argument
execute_called = False
my_exception_string = None

def add_subparser(subparser):
    parser = subparser.add_parser("mock_module", help="foo")
    parser.add_argument("output", help="foo")
    parser.add_argument("--force", action='store_true', help="foo")

def execute(args, execution_context):
    global execute_called
    execute_called = True
    if my_exception_string:
        raise Exception(my_exception_string)
