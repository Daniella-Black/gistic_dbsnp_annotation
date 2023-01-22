#!/usr/local/bin/python3

import pandas as pd
from dataclasses import dataclass

@dataclass
class SequenceRange:
    """Class for the start and end of a range."""
    name: str
    start: int
    end: int
    chrom: str
    def overlaps(self, other: "SequenceRange") -> bool:
        if self.chrom != other.chrom:
            return False
        return (other.start <= self.start <= other.end) or (other.start <= self.end <= other.end) or (self.start <= other.start <= self.end) or (self.start <= other.end <= self.end)
        
        
with open('readme.txt', 'w') as f:
    f.write('readme')
