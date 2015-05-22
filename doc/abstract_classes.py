"""This file contains abstract classes and methods which reflect the standard
architecture of Jacquard.

The following signatures are outlined:
    * FORMAT tags
    * variant callers
    * commands
"""


#pylint: disable=pointless-string-statement
import abc


"""Abstract class outlining requirements for adding a new tag to Jacquard"""
class NewTag(object):
    #pylint:disable=too-few-public-methods,abstract-class-not-used
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def add_tag_values(self, vcf_record):
        """Adds format tag and sample values to a VCF record"""



"""Abstract class outlining requirements for adding a new variant caller to
   Jacquard. The claim() method in this class calls _NewVcfReader()"""
class NewVariantCaller(object):
    #pylint:disable=too-few-public-methods,abstract-class-not-used
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def claim(self,file_readers):
        """Recognizes and claims MuTect VCFs form the set of all input VCFs.

        Each defined caller has a chance to evaluate and claim all the incoming
        files as something that it can process.

        Args:
            file_readers: the collection of currently unclaimed files

        Returns:
            A tuple of unclaimed readers and NewVariantCallerReaders.
        """



"""Abstract class outlining requirements for adding a new variant
  caller-specific VcfReader object to Jacquard. This class is called by
  the claim() method in NewVariantCaller()."""
class _NewVcfReader(object):
    #pylint:disable=abstract-class-not-used
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    @property
    def caller_name(self):
        """returns the caller name"""

    @abc.abstractmethod
    @property
    def file_name(self):
        """returns the file name"""

    @abc.abstractmethod
    def open(self):
        """opens the vcf reader"""

    @abc.abstractmethod
    def close(self):
        """closes the vcf reader"""

    @abc.abstractmethod
    @staticmethod
    def expected_file_format():
        """returns a list of file groupings (eg ['snp', 'indel']) to be used
        to later validate that the files are paired together correctly"""

    @abc.abstractmethod
    @property
    def metaheaders(self):
        """adds and returns new metaheaders"""

    @abc.abstractmethod
    @property
    def column_header(self):
        """returns column header"""

    @abc.abstractmethod
    def tagged_vcf_records(self):
        """yields tagged vcf_records"""

    @abc.abstractmethod
    def vcf_records(self):
        """yields vcf_records"""



"""Abstract methods outlining requirements for adding a new command to
   Jacquard."""
#pylint: disable=unused-argument
@abc.abstractmethod
def validate_args(args):
    """validates command-line arguments"""

@abc.abstractmethod
def get_required_input_output_types():
    """Returns a tuple of input and output types.

    Possible values are:
         'file'
         'directory'
    """

@abc.abstractmethod
def report_prediction(args):
    """returns a list of output file names be be created"""

@abc.abstractmethod
def add_subparser(subparser):
    """uses argparse to add a sub-parser"""

@abc.abstractmethod
def execute(args, execution_context):
    """executes private methods"""
