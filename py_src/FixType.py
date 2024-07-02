from enum import Enum


class FixType(Enum):
    NoFixNeeded = (False, True, 0)
    FixAdd = (True, True, 2)
    FixTrim = (True, True, 1)
    FixReplace = (True, True, 3)
    FailedReplace = (True, False, 5)
    FailedBadSegment = (False, False, 4)
    FailedNoAlignment = (True, False, 6)

    def __new__(cls, fix_attempted, good, rank):
        obj = object.__new__(cls)
        obj._value_ = (fix_attempted, good, rank)
        return obj

    def __init__(self, fix_attempted, good, rank):
        self.fix_attempted = fix_attempted
        self.good = good
        self.rank = rank
