class FixType():
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
