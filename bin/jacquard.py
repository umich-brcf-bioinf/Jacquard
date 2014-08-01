#!/usr/bin/python2.7
import argparse
import os
import sys 

import pivot_variants
import rollup_genes
import tag 
import normalize_varscan
import filter_somatic
import merge
import jacquard_utils

def main(modules, arguments):
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, 
    description='''type 'Jacquard -h <subcommand>' for help on a specific subcommand''', 
    epilog="authors: Jessica Bene, Chris Gates 07/2014")
    parser.add_argument("-v", "--version", action='version', version="Jacquard v{0}\nSupported variant callers:\n\t{1}".format(jacquard_utils.__version__, "\n\t".join([key + " " + value for key, value in jacquard_utils.caller_versions.items()])))
    subparsers = parser.add_subparsers(title="subcommands", dest="subparser_name")
    
    module_dispatch = {}
    for module in modules:
        module.add_subparser(subparsers)
        module_dispatch[module.__name__] = module

    execution_context = [
                         "##jacquard.version={0}".format(jacquard_utils.__version__), 
                         "##jacquard.command={0}".format(" ".join(arguments)), 
                         "##jacquard.cwd={0}".format(os.path.dirname(os.getcwd()))]

    args = parser.parse_args(arguments)
    module_dispatch[args.subparser_name].execute(args, execution_context)


if __name__ == "__main__":
    main([pivot_variants, rollup_genes, tag, normalize_varscan, filter_somatic, merge], sys.argv[1:])

