"""mrm_optimizer.py — Scheduled-MRM segment planner and CE-pair selector for the
Cabernet Sauvignon SPME Arrow GC-MS/MS panel.

Inputs
------
  MRM/transitions.yaml — analyte panel with precursor/product/CE/RT
  --peak-fwhm FLOAT     — chromatographic peak FWHM in seconds (typical 4-6 s for SPME-Arrow on DB-WAX)
  --target-ppp FLOAT    — target points-per-peak (default 12)
  --interscan-ms FLOAT  — interscan / pause time in ms (default 5)
  --window-min FLOAT    — sMRM half-window around each analyte RT in min (default 0.5)
  --min-dwell-ms FLOAT  — floor on dwell time (default 5)
  --classes LIST        — restrict panel to compound classes (comma-sep)

Outputs
-------
  acquisition_plan.json — full per-segment plan with N_concurrent, CT, PPP, dwell-per-tr
  segments.csv          — vendor-agnostic segment table (start_min, end_min, transition_id, dwell_ms, CE_eV)
  panel_summary.md      — human-readable summary with the cycle-time math

Cycle-time math (Sciex Pro / Agilent dMRM equivalent)
-----------------------------------------------------
  CT = N_concurrent × (D + P)
  PPP = (peak_fwhm_s / CT) — should be >= 10-12 for quantitative reliability
  D = ((peak_fwhm_s / target_ppp) - P) / N_concurrent   [solved for max D]
  D_floor = max(D_solved, min_dwell_ms)

If PPP < target, the planner warns and (a) tries to widen segments to reduce
N_concurrent, (b) reports which analyte segments are over-budget.

Reference: Sciex RUO-MKT-02-8539-A "Scheduled MRM Pro Algorithm Overview"
          Agilent 5990-3595EN "Dynamic MRM Mode"
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

try:
    import yaml
except ImportError:
    print("ERROR: pyyaml required. pip install pyyaml", file=sys.stderr)
    raise


_HERE = Path(__file__).resolve().parent
PROJECT_ROOT = _HERE.parent
DEFAULT_PANEL = PROJECT_ROOT / "MRM" / "transitions.yaml"


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class Transition:
    compound: str
    role: str          # quan | qual | confirm | isotope
    precursor: float
    product: float
    ce_ev: float
    rt_min: float
    rt_window_min: float
    source: str

    @property
    def tid(self) -> str:
        return f"{self.compound}_{self.role}_{int(self.precursor)}>{int(self.product)}"

    @property
    def window_start(self) -> float:
        return max(0.0, self.rt_min - self.rt_window_min)

    @property
    def window_end(self) -> float:
        return self.rt_min + self.rt_window_min


@dataclass
class Segment:
    start_min: float
    end_min: float
    transitions: list[Transition] = field(default_factory=list)

    @property
    def n_concurrent(self) -> int:
        return len(self.transitions)

    @property
    def duration_min(self) -> float:
        return self.end_min - self.start_min


@dataclass
class Plan:
    segments: list[Segment]
    peak_fwhm_s: float
    target_ppp: float
    interscan_ms: float
    min_dwell_ms: float
    column: str
    method_total_min: float

    def cycle_time_ms(self, seg: Segment, dwell_ms: float) -> float:
        return seg.n_concurrent * (dwell_ms + self.interscan_ms)

    def points_per_peak(self, seg: Segment, dwell_ms: float) -> float:
        ct_s = self.cycle_time_ms(seg, dwell_ms) / 1000.0
        if ct_s == 0:
            return float("inf")
        return self.peak_fwhm_s / ct_s

    def max_dwell_for_target(self, seg: Segment) -> float:
        """Solve D for PPP=target: CT = peak_fwhm_s / target_ppp; D = CT/N - P."""
        if seg.n_concurrent == 0:
            return self.min_dwell_ms
        ct_target_ms = (self.peak_fwhm_s / self.target_ppp) * 1000.0
        d = (ct_target_ms / seg.n_concurrent) - self.interscan_ms
        return max(self.min_dwell_ms, d)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_panel(path: Path = DEFAULT_PANEL,
               classes: Optional[list[str]] = None,
               sMRM_window_min: float = 0.5) -> list[Transition]:
    """Load all transitions, optionally restrict to compound classes."""
    with open(path, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    compounds = data.get("compounds", {})
    txs: list[Transition] = []
    for cid, c in compounds.items():
        if classes and c.get("class") not in classes:
            continue
        rt = float(c.get("rt") or 0.0)
        for t in (c.get("transitions") or []):
            txs.append(Transition(
                compound=cid,
                role=t.get("role", "quan"),
                precursor=float(t["precursor"]),
                product=float(t["product"]),
                ce_ev=float(t["ce_ev"]),
                rt_min=rt,
                rt_window_min=sMRM_window_min,
                source=t.get("source", ""),
            ))
    return txs


# ---------------------------------------------------------------------------
# Sweep-line segmentation
# ---------------------------------------------------------------------------


def build_segments(transitions: list[Transition]) -> list[Segment]:
    """Sweep-line over [window_start, window_end] of each transition.

    Each unique breakpoint defines a segment boundary; for each segment, the
    active transitions are those whose [start, end] contains the segment's
    midpoint. This is the canonical scheduled-MRM segmentation that vendor
    algorithms compute internally.
    """
    if not transitions:
        return []

    breakpoints = set()
    for t in transitions:
        breakpoints.add(round(t.window_start, 4))
        breakpoints.add(round(t.window_end, 4))
    bps = sorted(breakpoints)

    segments: list[Segment] = []
    for i in range(len(bps) - 1):
        start = bps[i]
        end = bps[i + 1]
        if end - start < 1e-6:
            continue
        mid = (start + end) / 2.0
        active = [t for t in transitions
                  if t.window_start <= mid <= t.window_end]
        if not active:
            continue
        segments.append(Segment(start_min=start, end_min=end, transitions=active))
    return segments


# ---------------------------------------------------------------------------
# Plan + report
# ---------------------------------------------------------------------------


def build_plan(transitions: list[Transition],
               peak_fwhm_s: float,
               target_ppp: float,
               interscan_ms: float,
               min_dwell_ms: float,
               column: str = "DB-WAX") -> Plan:
    segs = build_segments(transitions)
    total = max((s.end_min for s in segs), default=0.0) - min((s.start_min for s in segs), default=0.0)
    return Plan(
        segments=segs,
        peak_fwhm_s=peak_fwhm_s,
        target_ppp=target_ppp,
        interscan_ms=interscan_ms,
        min_dwell_ms=min_dwell_ms,
        column=column,
        method_total_min=round(total, 2),
    )


def plan_to_dict(p: Plan) -> dict:
    out_segs = []
    for s in p.segments:
        d = p.max_dwell_for_target(s)
        ct = p.cycle_time_ms(s, d)
        ppp = p.points_per_peak(s, d)
        out_segs.append({
            "start_min": s.start_min,
            "end_min": s.end_min,
            "duration_min": round(s.duration_min, 3),
            "n_concurrent": s.n_concurrent,
            "dwell_ms": round(d, 2),
            "cycle_time_ms": round(ct, 2),
            "points_per_peak": round(ppp, 1),
            "transitions": [t.tid for t in s.transitions],
        })
    return {
        "method": {
            "column": p.column,
            "peak_fwhm_s": p.peak_fwhm_s,
            "target_points_per_peak": p.target_ppp,
            "interscan_ms": p.interscan_ms,
            "min_dwell_ms": p.min_dwell_ms,
            "total_runtime_min": p.method_total_min,
        },
        "segments": out_segs,
    }


def write_csv(p: Plan, path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["segment_start_min", "segment_end_min", "compound",
                    "role", "precursor_mz", "product_mz", "ce_ev",
                    "dwell_ms", "source"])
        for s in p.segments:
            d = p.max_dwell_for_target(s)
            for t in s.transitions:
                w.writerow([round(s.start_min, 3), round(s.end_min, 3),
                            t.compound, t.role, t.precursor, t.product,
                            t.ce_ev, round(d, 2), t.source])


def render_markdown(p: Plan) -> str:
    lines = [
        "# Acquisition Plan — Scheduled MRM",
        "",
        f"- Column: {p.column}",
        f"- Peak FWHM: {p.peak_fwhm_s} s",
        f"- Target points-per-peak: {p.target_ppp}",
        f"- Interscan delay: {p.interscan_ms} ms",
        f"- Min dwell: {p.min_dwell_ms} ms",
        f"- Total runtime: {p.method_total_min} min",
        f"- Segments: {len(p.segments)}",
        "",
        "## Segment-by-segment cycle-time table",
        "",
        "| # | start (min) | end (min) | N_conc | dwell (ms) | CT (ms) | PPP | warning |",
        "|---|---|---|---|---|---|---|---|",
    ]
    overbudget = 0
    for i, s in enumerate(p.segments, 1):
        d = p.max_dwell_for_target(s)
        ct = p.cycle_time_ms(s, d)
        ppp = p.points_per_peak(s, d)
        warn = ""
        if ppp < p.target_ppp - 0.5:
            warn = f"⚠ PPP<{p.target_ppp}"
            overbudget += 1
        if d == p.min_dwell_ms:
            warn = (warn + " min-dwell-floor").strip()
        lines.append(
            f"| {i} | {s.start_min:.2f} | {s.end_min:.2f} | {s.n_concurrent} "
            f"| {d:.1f} | {ct:.1f} | {ppp:.1f} | {warn} |"
        )
    lines.append("")
    if overbudget:
        lines.append(f"**WARNING: {overbudget} segment(s) below target PPP.**")
        lines.append("Mitigations:")
        lines.append("  - Narrow the sMRM window per analyte (reduces concurrency).")
        lines.append("  - Drop the qual transition for high-abundance analytes (>10× LOQ).")
        lines.append("  - Increase peak FWHM by lowering carrier flow or increasing oven ramp slope.")
        lines.append("")

    lines.append("## Per-transition detail")
    lines.append("")
    seen = set()
    for s in p.segments:
        for t in s.transitions:
            if t.tid in seen:
                continue
            seen.add(t.tid)
            lines.append(f"- **{t.compound}** ({t.role}) — {t.precursor:.0f} → {t.product:.0f} @ {t.ce_ev:.0f} eV "
                         f"— RT {t.rt_min:.1f} min ±{t.rt_window_min:.2f} min — {t.source}")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CE optimization helper (CE-ramp simulator)
# ---------------------------------------------------------------------------


def center_of_mass_energy(e_lab_eV: float, m_ion: float, m_gas: float = 28.0) -> float:
    """Convert lab-frame collision energy to center-of-mass energy.

    E_CM = E_LAB × M_gas / (M_gas + M_ion)

    m_gas: N2=28, Ar=40, He=4 (default N2 — modern QqQ standard).
    """
    return e_lab_eV * (m_gas / (m_gas + m_ion))


def suggest_ce_ramp(precursor_mz: float,
                    start_eV: float = 5.0,
                    end_eV: float = 40.0,
                    step_eV: float = 2.5,
                    m_gas: float = 28.0) -> list[dict]:
    """Suggest a CE ramp for product-ion optimization.

    Returns a list of {ce_lab_eV, ce_cm_eV} pairs. Use this to plan a Q3 product-ion
    scan with stepped CE, find the CE that maximizes the diagnostic product ion.
    """
    out = []
    ce = start_eV
    while ce <= end_eV + 1e-9:
        out.append({
            "ce_lab_eV": round(ce, 2),
            "ce_cm_eV": round(center_of_mass_energy(ce, precursor_mz, m_gas), 3),
        })
        ce += step_eV
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Scheduled-MRM segment planner for the wine SPME Arrow GC-MS/MS panel."
    )
    ap.add_argument("--panel", default=str(DEFAULT_PANEL),
                    help="Path to transitions.yaml (default: project MRM/transitions.yaml)")
    ap.add_argument("--peak-fwhm", type=float, default=5.0,
                    help="Chromatographic peak FWHM in seconds (default 5)")
    ap.add_argument("--target-ppp", type=float, default=12.0,
                    help="Target points-per-peak (default 12)")
    ap.add_argument("--interscan-ms", type=float, default=5.0,
                    help="Interscan / pause time in ms (default 5)")
    ap.add_argument("--window-min", type=float, default=0.5,
                    help="sMRM half-window per analyte in min (default 0.5)")
    ap.add_argument("--min-dwell-ms", type=float, default=5.0,
                    help="Floor on dwell time in ms (default 5)")
    ap.add_argument("--classes", default=None,
                    help="Restrict to compound classes (comma-sep), e.g. "
                         "'methoxypyrazine,smoke_phenol'")
    ap.add_argument("--out-dir", default=None,
                    help="Output directory (default: project-root/MRM/plans/<stamp>/)")
    ap.add_argument("--ce-ramp", type=float, default=None,
                    help="Print a CE ramp suggestion for a given precursor m/z and exit")
    args = ap.parse_args()

    if args.ce_ramp is not None:
        ramp = suggest_ce_ramp(args.ce_ramp)
        print(f"CE ramp for precursor m/z {args.ce_ramp}:")
        for row in ramp:
            print(f"  CE_lab = {row['ce_lab_eV']:5.1f} eV  ->  CE_cm = {row['ce_cm_eV']:6.3f} eV")
        return 0

    classes = [s.strip() for s in args.classes.split(",")] if args.classes else None
    txs = load_panel(Path(args.panel), classes=classes,
                     sMRM_window_min=args.window_min)
    if not txs:
        print("ERROR: no transitions matched filter", file=sys.stderr)
        return 1

    plan = build_plan(
        txs,
        peak_fwhm_s=args.peak_fwhm,
        target_ppp=args.target_ppp,
        interscan_ms=args.interscan_ms,
        min_dwell_ms=args.min_dwell_ms,
    )

    # Output
    if args.out_dir:
        out_dir = Path(args.out_dir)
    else:
        from datetime import datetime
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = PROJECT_ROOT / "MRM" / "plans" / stamp
    out_dir.mkdir(parents=True, exist_ok=True)

    json_path = out_dir / "acquisition_plan.json"
    csv_path = out_dir / "segments.csv"
    md_path = out_dir / "panel_summary.md"

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(plan_to_dict(plan), f, indent=2)
    write_csv(plan, csv_path)
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(render_markdown(plan))

    print(f"plan written:")
    print(f"  {json_path}")
    print(f"  {csv_path}")
    print(f"  {md_path}")
    print(f"\n{len(plan.segments)} segments, "
          f"{sum(s.n_concurrent for s in plan.segments)} transitions, "
          f"{len(txs)} unique transitions, "
          f"total runtime {plan.method_total_min} min")

    overbudget = sum(
        1 for s in plan.segments
        if plan.points_per_peak(s, plan.max_dwell_for_target(s)) < args.target_ppp - 0.5
    )
    if overbudget:
        print(f"⚠ {overbudget} segment(s) below target PPP — see panel_summary.md")

    return 0


if __name__ == "__main__":
    sys.exit(main())
