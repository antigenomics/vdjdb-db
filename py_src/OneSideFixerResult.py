from dataclasses import dataclass
from FixType import FixType


@dataclass
class OneSideFixerResult:
    cdr3: str
    segmentId: str
    FixType: FixType
    x: int = -1
