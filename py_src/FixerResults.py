from dataclasses import dataclass
from enum import Enum
from FixType import FixType


@dataclass
class FixerResult:
    cdr3: str
    cdr3_old: str
    vId: str
    jId: str
    vEnd: int
    jStart: int
    vCanonical: bool
    jCanonical: bool
    fixNeeded: bool
    vFixType: FixType
    jFixType: FixType

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
