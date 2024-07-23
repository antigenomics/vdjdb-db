class FixType:
    def __init__(self, fix_attempted, good, rank, name):
        self.fix_attempted = fix_attempted
        self.good = good
        self.rank = rank
        self.name = name


NoFixNeeded = FixType(False, True, 0, 'NoFixNeeded')
FixAdd = FixType(True, True, 2, 'FixAdd')
FixTrim = FixType(True, True, 1, 'FixTrim')
FixReplace = FixType(True, True, 3, 'FixReplace')
FailedReplace = FixType(True, False, 5, 'FailedReplace')
FailedBadSegment = FixType(False, False, 4, 'FailedBadSegment')
FailedNoAlignment = FixType(True, False, 6, 'FailedNoAlignment')
