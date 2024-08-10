from dataclasses import dataclass


@dataclass
class FixType:
    """
    class for marking different scenarios of fix applied
    """
    fix_attempted: bool
    good: bool
    rank: int
    name: str


NoFixNeeded = FixType(False, True, 0, 'NoFixNeeded')
FixAdd = FixType(True, True, 2, 'FixAdd')
FixTrim = FixType(True, True, 1, 'FixTrim')
FixReplace = FixType(True, True, 3, 'FixReplace')
FailedReplace = FixType(True, False, 5, 'FailedReplace')
FailedBadSegment = FixType(False, False, 4, 'FailedBadSegment')
FailedNoAlignment = FixType(True, False, 6, 'FailedNoAlignment')


@dataclass
class OneSideFixerResult:
    """
    dataclass for results of one side fix
    """
    cdr3: str
    segmentId: str
    FixType: FixType
    x: int = -1


class FixerResult:
    """
    dataclass for fixer results of whole CDR3
    """
    def __init__(self, cdr3: str, cdr3_old: str, fix_needed: bool, v_rnd: int, j_dtart: int,
                 v_id: str, v_fix_type: FixType, j_id: str, j_fix_type: FixType):
        self.cdr3 = cdr3
        self.cdr3_old = cdr3_old
        self.fixNeeded = fix_needed
        self.vId = v_id
        self.vEnd = v_rnd
        self.jStart = j_dtart
        self.vFixType = v_fix_type
        self.jFixType = j_fix_type
        self.jId = j_id
        self.vCanonical = cdr3.startswith("C")
        self.jCanonical = cdr3.endswith("F") or cdr3.endswith("W")

    def is_good(self):
        return self.vFixType.good and self.jFixType.good

    def results_to_dict(self):
        return {
            "cdr3": self.cdr3,
            "cdr3_old": self.cdr3_old,
            "fixNeeded": self.fixNeeded,
            "good": self.is_good(),
            "jCanonical": self.jCanonical,
            "jFixType": self.jFixType.name,
            "jId": self.jId,
            "jStart": self.jStart,
            "vCanonical": self.vCanonical,
            "vEnd": self.vEnd,
            "vFixType": self.vFixType.name,
            "vId": self.vId
        }
