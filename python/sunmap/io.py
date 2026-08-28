"""Input reading for the SUNMAP port.

sunmap.f reads map_data.PRN with FORMAT(5(I8)): five fixed-width
8-character integer fields -- ICSG, iday, ICWSA, ILNS, ICLD -- per
data record, preceded by a FORMAT(I8,I8,A75) batch header line
(record count, N, label). map_data.PRN itself isn't in this repo.

The closest available data is rgo_data_0508.prn / _0808.prn / _5809.prn
(derived from rgo_data.sql): same five leading columns plus a trailing
CUA column that sunmap.f ignores and recomputes itself (see model.py).
Their first line is a human-readable column-header row, not a numeric
NREC/N/LABEL batch header, so this reader treats the whole file as one
batch (skip the header row, parse every remaining row as a record) --
see model.py's docstring for why the interactive batching loop was
dropped.
"""

from __future__ import annotations

from pathlib import Path

from .model import SunspotRecord


def read_prn(path: str | Path) -> list[SunspotRecord]:
    """Read a rgo_data_*.prn file into a list of SunspotRecord.

    Skips a non-numeric first line (column-header row); parses ID, day,
    CWSA, LNS, CLD from each remaining line (a trailing CUA column, if
    present, is ignored -- sunmap.f recomputes CUA itself).
    """
    records = []
    with open(path) as f:
        lines = f.readlines()

    for line in lines:
        parts = line.replace(",", " ").split()
        if not parts:
            continue
        try:
            values = [int(float(v)) for v in parts[:5]]
        except ValueError:
            continue  # header row or other non-numeric line
        if len(values) < 5:
            continue
        csg, day, cwsa, lns, cld = values
        records.append(SunspotRecord(csg=csg, day=day, cwsa=cwsa, lns=lns, cld=cld))

    return records
