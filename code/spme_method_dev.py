"""spme_method_dev.py — SPME Arrow method-development toolkit.

Functions
---------
  plackett_burman(n_factors, levels=2)   — N-run PB screening matrix (N = next multiple of 4)
  box_behnken(n_factors)                 — Box-Behnken response-surface design
  central_composite(n_factors, alpha)    — face-centered or rotatable CCD
  fit_quadratic_rsm(X, y)                — full quadratic RSM regression with ANOVA
  desirability(values, lo, hi, weights)  — Derringer-Suich one-sided + two-sided composite
  pawliszyn_isotherm(t, t_eq, n_inf)     — first-order extraction-isotherm model
  fit_isotherm(times, signals)           — least-squares fit of pawliszyn_isotherm
  salting_out_constant(Cs, log_K_ratios) — Setschenow-equation fit

CLI
---
  python spme_method_dev.py screen     --factors 7
  python spme_method_dev.py optimize   --design bbd --factors 3
  python spme_method_dev.py isotherm   --data times_signals.csv
  python spme_method_dev.py desire     --responses 0.5,1.0,0.8 --los 0,0,0 --his 1,1,1

The CLI writes design matrices to MRM/plans/<stamp>/design_*.csv and fit
summaries to JSON.

Reference: Derringer & Suich 1980 J Qual Tech 12, 214; Pawliszyn 2012
SPME book; Welke 2024 OENO; Saha 2018 Foods 7:127.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from itertools import combinations, product
from pathlib import Path
from typing import Optional

_HERE = Path(__file__).resolve().parent
PROJECT_ROOT = _HERE.parent


# ---------------------------------------------------------------------------
# Plackett-Burman screening
# ---------------------------------------------------------------------------


# Hadamard-derived PB matrices for the most common sizes. Each row is the
# generating row of the cyclic Hadamard matrix; rows after are cyclic shifts,
# and the final row is all -1. These are the canonical 1946 Plackett-Burman
# matrices. Higher orders (28, 32, ...) follow Paley construction.
_PB_GENERATORS = {
    4:  "++-",
    8:  "+++-+--",
    12: "++-+++---+-",
    16: "++++-+-++--+---",
    20: "++--++++-+-+----++-",
    24: "+++++-+-++--++--+-+----",
}


def plackett_burman(n_factors: int) -> list[list[int]]:
    """Return an N-run PB matrix as list-of-rows; entries are ±1.

    N = smallest multiple of 4 with N > n_factors. Returns the first n_factors
    columns. The (N-1)-th column is the all-equal "fold-over" check column.
    """
    n_runs = ((n_factors // 4) + 1) * 4
    while n_runs not in _PB_GENERATORS and n_runs < 32:
        n_runs += 4
    if n_runs not in _PB_GENERATORS:
        raise ValueError(f"no PB matrix bundled for n_runs={n_runs}; "
                         f"use central_composite or supply Paley construction")
    gen = _PB_GENERATORS[n_runs]
    # The first row is the generator. Subsequent rows are right-cyclic shifts.
    # The last row is all -1.
    rows: list[list[int]] = []
    cur = [1 if c == "+" else -1 for c in gen]
    for _ in range(n_runs - 1):
        rows.append(cur[:])
        cur = [cur[-1]] + cur[:-1]
    rows.append([-1] * (n_runs - 1))
    # Keep only the first n_factors columns
    return [r[:n_factors] for r in rows]


# ---------------------------------------------------------------------------
# Box-Behnken design
# ---------------------------------------------------------------------------


def box_behnken(n_factors: int, n_center: int = 3) -> list[list[int]]:
    """Box-Behnken design for n_factors (typically 3-5), levels coded -1/0/+1.

    Construction: for each pair of factors (i, j), the 2² factorial points
    (±1, ±1) are run at those two factors with all other factors at 0. Plus
    n_center center points (all factors = 0). This is the canonical Box-Behnken
    construction (Box & Behnken 1960).
    """
    if n_factors < 3:
        raise ValueError("Box-Behnken requires at least 3 factors")
    rows: list[list[int]] = []
    for i, j in combinations(range(n_factors), 2):
        for a in (-1, 1):
            for b in (-1, 1):
                r = [0] * n_factors
                r[i] = a
                r[j] = b
                rows.append(r)
    for _ in range(n_center):
        rows.append([0] * n_factors)
    return rows


# ---------------------------------------------------------------------------
# Central composite design
# ---------------------------------------------------------------------------


def central_composite(n_factors: int, alpha: Optional[float] = None,
                      n_center: int = 3, face_centered: bool = False) -> list[list[float]]:
    """Central composite design.

    factorial: 2^n_factors (or fractional for n>=5)
    axial: 2*n_factors points at ±alpha on each axis with all others at 0
    center: n_center replicates at origin

    alpha:
       - rotatable: alpha = (n_factorial)^(1/4)  — typical
       - face-centered: alpha = 1                 — fits inside the cube
       - orthogonal: solved from formula          — not implemented; cite Myers-Montgomery
    """
    rows: list[list[float]] = []
    # Factorial portion (full 2^k)
    for combo in product((-1.0, 1.0), repeat=n_factors):
        rows.append(list(combo))
    if alpha is None:
        if face_centered:
            alpha = 1.0
        else:
            alpha = (2 ** n_factors) ** 0.25  # rotatable
    # Axial / star points
    for i in range(n_factors):
        for sign in (-alpha, alpha):
            r = [0.0] * n_factors
            r[i] = sign
            rows.append(r)
    # Center
    for _ in range(n_center):
        rows.append([0.0] * n_factors)
    return rows


# ---------------------------------------------------------------------------
# Quadratic RSM regression (no scipy/numpy dependency — Gauss elimination)
# ---------------------------------------------------------------------------


def _expand_quadratic(x: list[float]) -> list[float]:
    """Return the design row [1, x1, x2, ..., x1*x2, ..., x1², x2², ...]."""
    k = len(x)
    row = [1.0]
    row.extend(x)
    for i, j in combinations(range(k), 2):
        row.append(x[i] * x[j])
    for i in range(k):
        row.append(x[i] * x[i])
    return row


def _matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    m, k1 = len(a), len(a[0])
    k2, n = len(b), len(b[0])
    assert k1 == k2
    out = [[0.0] * n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            s = 0.0
            for kk in range(k1):
                s += a[i][kk] * b[kk][j]
            out[i][j] = s
    return out


def _transpose(a: list[list[float]]) -> list[list[float]]:
    return [list(r) for r in zip(*a)]


def _inv(a: list[list[float]]) -> list[list[float]]:
    """Matrix inverse via Gauss-Jordan elimination with partial pivoting."""
    n = len(a)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)]
           for i, row in enumerate(a)]
    for col in range(n):
        # Pivot
        pivot_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot_row][col]) < 1e-12:
            raise ValueError("singular matrix in _inv")
        aug[col], aug[pivot_row] = aug[pivot_row], aug[col]
        pv = aug[col][col]
        aug[col] = [v / pv for v in aug[col]]
        for r in range(n):
            if r == col:
                continue
            factor = aug[r][col]
            aug[r] = [v - factor * w for v, w in zip(aug[r], aug[col])]
    return [row[n:] for row in aug]


def fit_quadratic_rsm(X: list[list[float]], y: list[float]) -> dict:
    """Fit y = beta0 + sum(beta_i*x_i) + sum(beta_ij*x_i*x_j) + sum(beta_ii*x_i²).

    Returns {coefficients, predicted, residuals, R², adj_R², ANOVA p-value placeholder}.
    """
    n = len(y)
    Xd = [_expand_quadratic(row) for row in X]
    Xt = _transpose(Xd)
    XtX = _matmul(Xt, Xd)
    XtY = _matmul(Xt, [[v] for v in y])
    try:
        XtX_inv = _inv(XtX)
    except ValueError as e:
        return {"error": str(e)}
    beta = _matmul(XtX_inv, XtY)
    beta = [row[0] for row in beta]
    y_hat = []
    for row in Xd:
        y_hat.append(sum(b * r for b, r in zip(beta, row)))
    residuals = [yi - yh for yi, yh in zip(y, y_hat)]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    ss_res = sum(r * r for r in residuals)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    p = len(beta)
    adj_r_squared = 1.0 - (1.0 - r_squared) * (n - 1) / max(1, (n - p))
    return {
        "n": n,
        "n_terms": p,
        "coefficients": beta,
        "predicted": y_hat,
        "residuals": residuals,
        "r_squared": round(r_squared, 4),
        "adj_r_squared": round(adj_r_squared, 4),
        "rss": round(ss_res, 4),
    }


# ---------------------------------------------------------------------------
# Desirability function (Derringer & Suich 1980)
# ---------------------------------------------------------------------------


def desirability_max(value: float, lo: float, hi: float, r: float = 1.0) -> float:
    """One-sided desirability — larger is better.

    d = 0 if value <= lo
    d = ((value - lo) / (hi - lo)) ** r  if lo < value < hi
    d = 1 if value >= hi
    """
    if value <= lo:
        return 0.0
    if value >= hi:
        return 1.0
    return ((value - lo) / (hi - lo)) ** r


def desirability_min(value: float, lo: float, hi: float, r: float = 1.0) -> float:
    """One-sided desirability — smaller is better."""
    if value <= lo:
        return 1.0
    if value >= hi:
        return 0.0
    return ((hi - value) / (hi - lo)) ** r


def desirability_target(value: float, lo: float, target: float, hi: float,
                         r_low: float = 1.0, r_high: float = 1.0) -> float:
    if value < lo or value > hi:
        return 0.0
    if value <= target:
        return ((value - lo) / (target - lo)) ** r_low
    return ((hi - value) / (hi - target)) ** r_high


def composite_desirability(individual_d: list[float],
                            weights: Optional[list[float]] = None) -> float:
    """Geometric mean of individual desirabilities (weighted)."""
    if not individual_d:
        return 0.0
    if weights is None:
        weights = [1.0] * len(individual_d)
    w_sum = sum(weights)
    if w_sum == 0:
        return 0.0
    log_d = 0.0
    for d, w in zip(individual_d, weights):
        if d <= 0:
            return 0.0  # any zero kills the geometric mean
        log_d += w * math.log(d)
    return math.exp(log_d / w_sum)


# ---------------------------------------------------------------------------
# SPME extraction-isotherm model (Pawliszyn 1990, three-phase)
# ---------------------------------------------------------------------------


def pawliszyn_isotherm(t: float, t_eq: float, n_inf: float) -> float:
    """First-order extraction isotherm (Pawliszyn 1990 simplification).

    n(t) = n_inf * (1 - exp(-t / tau))
    where t_eq is the 95% equilibrium time, and tau = t_eq / 3.0
    (since 1 - exp(-3) ≈ 0.95).
    """
    tau = t_eq / 3.0
    return n_inf * (1.0 - math.exp(-t / tau))


def fit_isotherm(times: list[float], signals: list[float],
                 max_iter: int = 200) -> dict:
    """Levenberg-Marquardt-style fit of (t_eq, n_inf) via grid + golden-section.

    Returns {t_eq, n_inf, rmse, r_squared}.
    """
    if len(times) < 3:
        return {"error": "need >= 3 datapoints"}

    n_inf_init = max(signals) * 1.1
    t_eq_init = times[-1] * 0.8

    def loss(t_eq: float, n_inf: float) -> float:
        s = 0.0
        for t, y in zip(times, signals):
            yh = pawliszyn_isotherm(t, t_eq, n_inf)
            s += (y - yh) ** 2
        return s

    # Alternating-axis golden-section search
    t_eq = t_eq_init
    n_inf = n_inf_init
    for _ in range(max_iter):
        # 1D search over t_eq
        lo, hi = max(0.1, t_eq * 0.1), t_eq * 5.0
        for _ in range(40):
            phi = (math.sqrt(5) - 1) / 2
            x1 = hi - phi * (hi - lo)
            x2 = lo + phi * (hi - lo)
            if loss(x1, n_inf) < loss(x2, n_inf):
                hi = x2
            else:
                lo = x1
        new_t_eq = (lo + hi) / 2
        # 1D search over n_inf
        lo, hi = max(signals) * 0.5, max(signals) * 3.0
        for _ in range(40):
            phi = (math.sqrt(5) - 1) / 2
            x1 = hi - phi * (hi - lo)
            x2 = lo + phi * (hi - lo)
            if loss(new_t_eq, x1) < loss(new_t_eq, x2):
                hi = x2
            else:
                lo = x1
        new_n_inf = (lo + hi) / 2
        if abs(new_t_eq - t_eq) < 1e-6 and abs(new_n_inf - n_inf) < 1e-6:
            break
        t_eq, n_inf = new_t_eq, new_n_inf

    # Stats
    rss = loss(t_eq, n_inf)
    rmse = math.sqrt(rss / len(times))
    y_mean = sum(signals) / len(signals)
    ss_tot = sum((y - y_mean) ** 2 for y in signals)
    r_squared = 1.0 - rss / ss_tot if ss_tot > 0 else 1.0
    return {
        "t_eq_min": round(t_eq, 3),
        "n_inf": round(n_inf, 3),
        "tau_min": round(t_eq / 3, 3),
        "rmse": round(rmse, 4),
        "r_squared": round(r_squared, 4),
        "n_points": len(times),
    }


# ---------------------------------------------------------------------------
# Setschenow / salting-out fit
# ---------------------------------------------------------------------------


def salting_out_fit(C_salt_M: list[float], log_K_app: list[float]) -> dict:
    """Linear fit log K_app = log K0 + k_S * [salt].

    k_S is the Setschenow constant (per molar). Positive k_S → salting-out;
    negative → salting-in.
    """
    n = len(C_salt_M)
    if n != len(log_K_app) or n < 2:
        return {"error": "need matching arrays of length >= 2"}
    mx = sum(C_salt_M) / n
    my = sum(log_K_app) / n
    num = sum((x - mx) * (y - my) for x, y in zip(C_salt_M, log_K_app))
    den = sum((x - mx) ** 2 for x in C_salt_M)
    if den == 0:
        return {"error": "zero variance in C_salt"}
    k_S = num / den
    log_K0 = my - k_S * mx
    y_hat = [log_K0 + k_S * x for x in C_salt_M]
    ss_tot = sum((y - my) ** 2 for y in log_K_app)
    ss_res = sum((y - yh) ** 2 for y, yh in zip(log_K_app, y_hat))
    return {
        "log_K0": round(log_K0, 4),
        "k_S_per_M": round(k_S, 4),
        "r_squared": round(1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0, 4),
        "interpretation": "salting-out" if k_S > 0 else ("salting-in" if k_S < 0 else "neutral"),
    }


# ---------------------------------------------------------------------------
# Default DoE factor catalog (for the dissertation)
# ---------------------------------------------------------------------------


_DEFAULT_FACTORS = {
    "extraction_temperature_C": {"low": 30, "high": 60, "unit": "°C", "rationale": "60°C is hard ceiling for IBMP artifact (Antalick)"},
    "extraction_time_min":     {"low": 15, "high": 60, "unit": "min", "rationale": "covers pre-equilibrium to full-equilibrium for most volatiles"},
    "NaCl_mass_g":              {"low": 0.5, "high": 4.0, "unit": "g", "rationale": "0-300 g/L; salting-out effect saturates ~250 g/L"},
    "sample_volume_mL":         {"low": 5, "high": 10, "unit": "mL", "rationale": "headspace ratio in 20 mL vial"},
    "ethanol_dilution":         {"low": 0, "high": 5, "unit": "x dilution", "rationale": "matrix dilution reduces ethanol depressive effect"},
    "agitation_rpm":            {"low": 250, "high": 500, "unit": "rpm", "rationale": "CTC PAL3 nominal 250-500"},
    "incubation_time_min":      {"low": 5, "high": 15, "unit": "min", "rationale": "thermal equilibration before fiber exposure"},
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _write_csv(rows: list[list[float]], header: list[str], path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def cmd_screen(args) -> int:
    factor_keys = (args.factors_named or "extraction_temperature_C,extraction_time_min,"
                                          "NaCl_mass_g,sample_volume_mL,ethanol_dilution,"
                                          "agitation_rpm,incubation_time_min").split(",")
    n = len(factor_keys)
    rows = plackett_burman(n)
    # Decode -1/+1 to actual levels using _DEFAULT_FACTORS if known
    header = ["run"] + factor_keys
    out = []
    for i, r in enumerate(rows, 1):
        decoded: list[float] = []
        for k, code in zip(factor_keys, r):
            f = _DEFAULT_FACTORS.get(k, {"low": -1, "high": 1})
            v = f["low"] if code == -1 else f["high"]
            decoded.append(v)
        out.append([i] + decoded)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = PROJECT_ROOT / "MRM" / "plans" / stamp
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "design_plackett_burman.csv"
    _write_csv(out, header, csv_path)
    print(f"wrote Plackett-Burman design ({n} factors, {len(rows)} runs) -> {csv_path}")
    return 0


def cmd_optimize(args) -> int:
    if args.design == "bbd":
        if args.factors < 3:
            print("Box-Behnken needs >= 3 factors; falling back to CCD", file=sys.stderr)
            return cmd_optimize(argparse.Namespace(design="ccd",
                                                   factors=args.factors,
                                                   alpha=args.alpha))
        rows = box_behnken(args.factors)
    elif args.design == "ccd":
        rows = central_composite(args.factors, alpha=args.alpha,
                                  face_centered=(args.alpha == "face"))
    else:
        print(f"unknown design '{args.design}'", file=sys.stderr)
        return 2

    factor_keys = (args.factors_named or
                   ",".join(list(_DEFAULT_FACTORS.keys())[:args.factors])).split(",")
    header = ["run"] + [f"x{i+1}_{k}" for i, k in enumerate(factor_keys)]

    out = []
    for i, r in enumerate(rows, 1):
        decoded = []
        for k, code in zip(factor_keys, r):
            f = _DEFAULT_FACTORS.get(k, {"low": -1, "high": 1})
            mid = (f["low"] + f["high"]) / 2
            half = (f["high"] - f["low"]) / 2
            v = mid + code * half
            decoded.append(round(v, 3))
        out.append([i] + decoded)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = PROJECT_ROOT / "MRM" / "plans" / stamp
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"design_{args.design}.csv"
    _write_csv(out, header, csv_path)
    print(f"wrote {args.design} design ({args.factors} factors, {len(rows)} runs) -> {csv_path}")
    return 0


def cmd_isotherm(args) -> int:
    rows = []
    with open(args.data, encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            if len(row) >= 2:
                try:
                    rows.append((float(row[0]), float(row[1])))
                except ValueError:
                    continue
    times = [r[0] for r in rows]
    signals = [r[1] for r in rows]
    fit = fit_isotherm(times, signals)
    print(json.dumps(fit, indent=2))
    return 0


def cmd_desire(args) -> int:
    vals = [float(x) for x in args.responses.split(",")]
    los = [float(x) for x in args.los.split(",")]
    his = [float(x) for x in args.his.split(",")]
    if not (len(vals) == len(los) == len(his)):
        print("responses/los/his must have equal lengths", file=sys.stderr)
        return 2
    weights = ([float(x) for x in args.weights.split(",")]
               if args.weights else [1.0] * len(vals))
    direction = (args.directions.split(",") if args.directions
                 else ["max"] * len(vals))
    d_each = []
    for v, lo, hi, dr in zip(vals, los, his, direction):
        if dr == "max":
            d_each.append(desirability_max(v, lo, hi))
        elif dr == "min":
            d_each.append(desirability_min(v, lo, hi))
        else:
            print(f"unknown direction '{dr}'", file=sys.stderr)
            return 2
    composite = composite_desirability(d_each, weights)
    print(json.dumps({
        "individual_d": [round(d, 4) for d in d_each],
        "weights": weights,
        "composite_D": round(composite, 4),
    }, indent=2))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description="SPME Arrow method-development toolkit (DoE + RSM + isotherm)."
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_screen = sub.add_parser("screen", help="Plackett-Burman screening design")
    p_screen.add_argument("--factors", type=int, default=7)
    p_screen.add_argument("--factors-named", default=None)
    p_screen.set_defaults(func=cmd_screen)

    p_opt = sub.add_parser("optimize", help="Box-Behnken or CCD optimization design")
    p_opt.add_argument("--design", choices=["bbd", "ccd"], default="bbd")
    p_opt.add_argument("--factors", type=int, default=3)
    p_opt.add_argument("--alpha", default=None,
                        help="CCD alpha: 'face' (1.0), 'rotatable' (default), or float")
    p_opt.add_argument("--factors-named", default=None)
    p_opt.set_defaults(func=cmd_optimize)

    p_iso = sub.add_parser("isotherm", help="Fit Pawliszyn first-order isotherm to time/signal data")
    p_iso.add_argument("--data", required=True, help="CSV with columns time_min,signal")
    p_iso.set_defaults(func=cmd_isotherm)

    p_des = sub.add_parser("desire", help="Composite desirability score (Derringer-Suich)")
    p_des.add_argument("--responses", required=True,
                        help="Comma-separated response values, e.g. '1500,8000,2.3e6'")
    p_des.add_argument("--los", required=True, help="Lower bounds")
    p_des.add_argument("--his", required=True, help="Upper bounds")
    p_des.add_argument("--weights", default=None)
    p_des.add_argument("--directions", default=None,
                        help="Per-response 'max' or 'min' (default all max)")
    p_des.set_defaults(func=cmd_desire)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
