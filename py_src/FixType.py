from enum import Enum


class FixType(Enum):

    def __new__(cls, fix_attempted, good, rank):
        obj = object.__new__(cls)
        obj._value_ = (fix_attempted, good, rank)
        return obj

    def __init__(self, fix_attempted, good, rank):
        self.fix_attempted = fix_attempted
        self.good = good
        self.rank = rank


NoFixNeeded = FixType(False, True, 0)
FixAdd = FixType(True, True, 2)
FixTrim = FixType(True, True, 1)
FixReplace = FixType(True, True, 3)
FailedReplace = FixType(True, False, 5)
FailedBadSegment = FixType(False, False, 4)
FailedNoAlignment = FixType(True, False, 6)
