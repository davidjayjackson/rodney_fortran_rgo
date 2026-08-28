"""Python port of the main SUNMAP record-processing loop (sunmap.f).

Scope
-----
This port covers PROGRAM SUNMAP's record loop (sunmap.f lines 96-198)
and the EXPO evolution-model chain (expo.py) it calls. Deliberately
NOT ported, per project scope decision:

- SPUR / ZCELLS: a self-contained SEIR-style polygon simulation whose
  required state files (POLYPARA.DAT, TIMENEW.DAT, TIMEOLD.DAT) are not
  present in this repository.
- moves / ran2: the legacy PRNG that only feeds SPUR.
- The RNSTNR-driven "random data" mode (ADAT != 'Y' branch): RNSTNR is
  never defined anywhere in the original source tree, so there is
  nothing to port -- this port only supports the file-driven path.
- The interactive "Continue (1), Stop (0)" batching prompt: that's
  terminal-UI plumbing, not part of the computation. `process_batch`
  below processes one full batch per call; drive it from a loop, a
  Quarto parameter, or however many batches you have.

ZCELLS stub
-----------
Where the original calls ZCELLS to fill INF(1..8) with a "visible
pores" estimate per polygon, this port substitutes 8 zeros
(`_zcells_stub`). That makes POSQ == 0 for every record where its
position in the batch is <= K (the "duration of umbra in days"
parameter) -- exactly what you'd get running the original without real
POLYPARA.DAT/TIMENEW.DAT/TIMEOLD.DAT data. EXPO is still invoked at
that point for call-shape parity with the source (see expo.py's
docstring: its result isn't consumed further in that branch, in the
original either).
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field

from .expo import expo

# "nearest neighbor, use later for map clusters" (sunmap.f line 81)
Z = 8


@dataclass
class SunmapParams:
    """The six values sunmap.f prompts for interactively (lines 46-73).

    Field names match the Fortran variable names so the two can be read
    side by side; see sunmap.f's prompt text for what each represents.
    """

    beta: float   # "Enter number of Umbra seen"
    delta: float  # "Percentage of Umbra with Penumbra"
    rr: float     # "Enter number of total sunspots seen"
    g: float      # "Percentage of sunspots in Umbra - gamma"
    p: float      # "Estimate percent of Umbra in groups P(0<P<1.0)"
    k: int        # "Estimate duration of Umbra in days - K(1<K<85)"
    alpha: float = 0.95  # hardcoded in the original (sunmap.f line 78)

    pr: float = field(init=False)
    x: float = field(init=False)

    def __post_init__(self):
        # sunmap.f line 58: rr = ((rr + 1) / 100.)
        self.rr = (self.rr + 1) / 100.0
        # sunmap.f line 66: Pr = (1.-P)/10
        self.pr = (1.0 - self.p) / 10.0
        # sunmap.f line 80: x = 1./Pr ("global frequency")
        self.x = 1.0 / self.pr if self.pr != 0 else math.inf


@dataclass
class SunspotRecord:
    """One data row: ICSG, iday, ICWSA, ILNS, ICLD (sunmap.f FORMAT(5(I8)))."""

    csg: int   # Greenwich/NOAA sunspot group number
    day: int   # day of month
    cwsa: int  # corrected whole spot area
    lns: int   # latitude
    cld: int   # Carrington longitude


def _zcells_stub(record_index: int) -> list[float]:
    """Stand-in for CALL ZCELLS(IREC-1, INF). See module docstring."""
    return [0.0] * Z


def process_batch(records: list[SunspotRecord], params: SunmapParams) -> list[dict]:
    """Port of the DO 100 IREC = NREC, NBATCH loop body (sunmap.f 111-175).

    `records` is one batch (N in the original = len(records)).  Returns
    one dict per record with the same fields sunmap.f writes to
    SUN.CSV: day, SUSQI, evolution (POSQ*SUSQS), CUA, ICWSA.
    """
    n = len(records)
    if n == 0:
        return []

    beta, delta, rr, g, p, k, alpha = (
        params.beta, params.delta, params.rr, params.g,
        params.p, params.k, params.alpha,
    )
    pr, x = params.pr, params.x

    last_day = 1
    # REMQ is only reassigned on even IREC (sunmap.f lines 160-164) and
    # keeps its prior value on odd IREC -- a carry-over quirk of the
    # original, preserved here rather than "fixed".
    remq = 0.0
    rows = []

    for irec, rec in enumerate(records, start=1):
        product = float(rec.lns) * float(rec.cld)
        if product < 0:
            # LNS (latitude) is signed -- Fortran's DSQRT of a negative
            # number yields NaN rather than raising, so this port does too
            # (Python's math.sqrt raises ValueError on a negative input).
            warnings.warn(
                f"record {irec} (day={rec.day}, csg={rec.csg}): "
                "LNS*CLD < 0, sqrt() would be undefined -- treated as NaN"
            )
            cua = math.nan
        elif rec.cwsa == 0:
            warnings.warn(
                f"record {irec} (day={rec.day}, csg={rec.csg}): ICWSA=0, "
                "CUA would divide by zero -- treated as NaN"
            )
            cua = math.nan
        else:
            cua = math.sqrt(product) / rec.cwsa

        if rec.day > last_day:
            last_day = rec.day

        susq = cua
        rinf = rr * cua
        freq0 = g * cua
        c = beta * cua
        pg = p * cua
        if freq0 == 0:
            warnings.warn(f"record {irec}: freq (g*CUA) is zero -- treated as NaN")
            freq = math.nan
        else:
            freq = (c / freq0) ** 2

        time = (float(irec) / n) + 0.01

        if irec <= k:
            inf = _zcells_stub(irec - 1)
            f = sum(inf)
            posq = f / Z
            expo(n, time, pg, rinf, freq, c)  # computed for parity; unused, see docstring
        else:
            posq = 1.0 / float(irec)

        if (irec + 1) % 2 == 1:
            remq = -(1.0 / float(irec))
            remq = remq * (g * (1 - alpha))
            remq = remq * delta
        # else: remq carries over from the previous record, matching sunmap.f

        susqs_denom = Z + pr * rr * (x * susq + x * remq)
        susqi_denom = Z + pr * beta * x * posq
        susqs = ((1 - pr) * rr * (n * susq + n * remq)) / susqs_denom if susqs_denom != 0 else math.nan
        susqi = ((1 - pr) * beta * n * posq) / susqi_denom if susqi_denom != 0 else math.nan
        susqs *= delta
        susqi *= (delta + g * alpha)

        rows.append({
            "day": rec.day,
            "SUSQI": susqi,
            "evolution": posq * susqs,
            "CUA": cua,
            "ICWSA": rec.cwsa,
        })

    return rows
