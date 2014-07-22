#!/usr/bin/python2.7
__version__ = 0.1

import argparse
import os
import sys 

import pivot_variants
import rollup_genes
import tag 
import normalize_varscan
import filter_somatic

def main(modules, arguments):
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, 
    description='''type 'Jacquard -h <subcommand>' for help on a specific subcommand''', 
    epilog="authors: Jessica Bene, Chris Gates 07/2014")
    parser.add_argument("-v", "--version", action='version', version="v{0}".format(__version__))
    subparsers = parser.add_subparsers(title="subcommands", dest="subparser_name")
    
    module_dispatch = {}
    for module in modules:
        module.add_subparser(subparsers)
        module_dispatch[module.__name__] = module

    execution_context = [
                         "##jacquard.version={0}".format(__version__), 
                         "##jacquard.command={0}".format(" ".join(arguments)), 
                         "##jacquard.cwd={0}".format(os.path.dirname(os.getcwd()))]

    args = parser.parse_args(arguments)
    module_dispatch[args.subparser_name].execute(args, execution_context)


if __name__ == "__main__":
    main([pivot_variants, rollup_genes, tag, normalize_varscan, filter_somatic], sys.argv[1:])

