from jacquard.variant_callers.varscan import Varscan
from jacquard.variant_callers.strelka import Strelka
from jacquard.variant_callers.mutect import Mutect
from jacquard.variant_callers.unknown import Unknown

class Factory():
    def __init__(self):
        self._callers = [Varscan(), Strelka(), Mutect(), Unknown()]
        
    def get_caller(self):