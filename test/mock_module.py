# pylint: disable=C0103,C0301,R0903,R0904,W0613,W0603
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
