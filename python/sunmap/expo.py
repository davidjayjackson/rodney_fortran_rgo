"""Python port of the EXPO evolution-model chain from sunmap.f.

Covers SUBROUTINE EXPO and its helpers ENERGY, MASS, JTEMP, E2, CPI.
This is a faithful, function-for-function port -- variable names and
the formula shapes match the Fortran source (sunmap.f lines 203-437)
so the two can be read side by side.

Not ported (out of scope, see python/sunmap/model.py docstring):
SPUR, ZCELLS, moves, ran2, and the RNSTNR call the original made in
its "random data" mode (RNSTNR is not defined anywhere in the
original source tree, so there is nothing to port).
"""

from __future__ import annotations

import math
import warnings

# Fortran source uses this truncated constant (PI = 3.141592654) rather
# than a full double-precision pi -- kept exactly for numerical parity.
PI = 3.141592654


def _safe_div(numerator: float, denominator: float, where: str) -> float:
    """Division with a warning instead of a crash on a zero denominator.

    The original Fortran has no such guard -- a zero denominator there
    would produce a runtime floating-point exception or an Infinity,
    depending on compiler/flags. Returning NaN here lets a batch keep
    processing instead of aborting, while still making the bad record
    visible (via the warning, and via NaN propagating into the output).
    """
    if denominator == 0:
        warnings.warn(f"{where}: division by zero, result set to NaN")
        return math.nan
    return numerator / denominator


def energy(time: float) -> float:
    """Port of SUBROUTINE ENERGY (sunmap.f lines 289-315).

    Original formula comment:
        h = exp((((ln(v)) * pi)-1)cos-1)
        h = 1/(pi*2)*h      <- present in the comment but commented out
                                in the actual code (line 312); NOT applied
                                here either, to match the executed Fortran.
    """
    lnh = math.log(time)
    cosrad = _safe_div(1.0, lnh * PI, "energy: cosrad")
    sinrad = abs(cosrad**2 - 1.0)
    temp = math.sqrt(sinrad)
    tanrad = _safe_div(temp, cosrad, "energy: tanrad")
    angle = -math.atan(tanrad)
    degtan = -angle
    exfunc = math.exp(degtan)
    h = _safe_div(1.0, exfunc, "energy: h")
    return h


def mass(time: float) -> float:
    """Port of SUBROUTINE MASS (sunmap.f lines 327-346)."""
    lnm = math.log(time)
    cosdeg = _safe_div(1.0, lnm / 0.4, "mass: cosdeg")
    tan_ = _safe_div(math.sqrt(abs(1 - cosdeg**2)), cosdeg, "mass: tan")
    arctan = math.atan(tan_)
    degree = arctan
    m = math.exp(degree)
    return m


def jtemp(time: float) -> float:
    """Port of SUBROUTINE JTEMP (sunmap.f lines 355-378).

    Guards against log(time) <= 0 ("FUNCTION FAILS AT ABSOLUTE 0" in the
    original comment) by returning 0.0, matching the GOTO 100 short
    circuit in the Fortran.
    """
    lnj = math.log(time)
    if lnj <= 0:
        return 0.0
    tanrad = math.exp(PI * math.log(lnj))
    angle = math.atan(tanrad)
    degtan = angle
    exfunc = math.exp(degtan)
    return exfunc


def e2(time: float) -> float:
    """Port of SUBROUTINE E2 (sunmap.f lines 385-417)."""
    cosrad = _safe_div(1.0, math.log(time) / 0.4, "e2: cosrad base") ** 2
    sinrad = abs(cosrad**2 - 1.0)
    temp = math.sqrt(sinrad)
    tanrad = _safe_div(temp, cosrad, "e2: tanrad")
    gm = _safe_div(1.0, math.atan(tanrad), "e2: gm")
    lne2 = math.log(PI * gm)
    tanrad2 = _safe_div(1.0, lne2, "e2: tanrad2")
    angle = math.atan(tanrad2)
    degtan = angle
    exfunc = math.exp(degtan)
    esqr = _safe_div(1.0, exfunc, "e2: esqr")
    return esqr


def cpi() -> float:
    """Port of SUBROUTINE CPI (sunmap.f lines 426-437).

    Not called anywhere in the active SUNMAP path (the main loop derives
    `c` from beta*CUA / beta*SUSQ instead) -- ported for completeness
    since it's part of the EXPO formula family described in the header
    comment.
    """
    return math.exp(math.exp((1.0 / PI) * 10.0))


def expo(n: float, time: float, pg: float, rinf: float, freq: float, c: float) -> float:
    """Port of SUBROUTINE EXPO (sunmap.f lines 203-278).

    Returns F, the SICe-derived evolution estimate. `n` is unused in the
    formula body in the original too -- kept as a parameter only for
    call-site parity with the Fortran signature.

    Caller note: in sunmap.f's main loop, this function's return value
    is computed but not read again in the same iteration -- POSQ (not F)
    drives the rest of that iteration's calculation (see model.py). This
    port preserves that call shape; it's on the caller whether to use F.
    """
    # "assume a solar minimum, longer latency period, DE-NORMALIZE for
    #  solar maximum, short latency period" (sunmap.f line 260-262)
    pg = 1.0 - _safe_div(pg - rinf, pg + rinf, "expo: pg normalization")
    h = energy(time)
    m = mass(time)
    esqr = e2(time)
    j = jtemp(time + 1)
    k = j**3
    kpg = k * pg
    if kpg > 0:
        kt = math.log10(kpg)
        if kt <= 0.0:
            kt = 5 / 2
    else:
        # DLOG10 of a non-positive argument would be a domain error in
        # Fortran; the original's IF(KT.LE.0.0) guard runs too late to
        # catch that case, so this port heads it off here instead.
        kt = 5 / 2
    s1 = _safe_div(1.0, 3 * PI * (h**3 * freq**3), "expo: s1")
    # Fortran: S2 = (3*DSQRT(ESQR)/4*PI*M*freq)**KT
    # * and / share precedence and evaluate left-to-right in Fortran, so
    # this is ((3*sqrt(esqr))/4)*PI*M*freq, NOT 3*sqrt(esqr)/(4*PI*M*freq)
    # -- despite what the header-comment formula (item 6) suggests. Kept
    # exactly as the code executes, not as the comment describes it.
    s2_base = (3 * math.sqrt(esqr) / 4) * PI * m * freq
    if s2_base < 0:
        # Fortran real**real exponentiation is EXP(KT*LOG(base)) under the
        # hood -- a negative base with non-integer KT is a domain error
        # there too. Guarded to NaN here instead of raising/going complex.
        warnings.warn("expo: s2_base negative with fractional exponent, result set to NaN")
        s2 = math.nan
    else:
        s2 = s2_base**kt
    sice = s1 * s2
    f = (sice * freq) * (m * c)
    return f
