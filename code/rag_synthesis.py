"""rag_synthesis.py — run dissertation-section queries against pgvector farm.

Writes results to projects/wine-spme-dissertation/_research/rag_synthesis.md
in a structured format that downstream CC-PRO composition can consume.

Each query corresponds to a dissertation section. Top-k results land with
similarity score, title, year, and (if available) distill TLDR card.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

# Force project repo dir
LIFE_ROOT = Path(r"C:/Users/jtole/Documents/2026/life")
os.chdir(LIFE_ROOT)
sys.path.insert(0, str(LIFE_ROOT))

from dotenv import load_dotenv
load_dotenv(LIFE_ROOT / ".env")

from research.embed.search import search


# (section_heading, query, k)
QUERIES = [
    ("2. Introduction — Cab Sauv aroma complexity",
     "Cabernet Sauvignon volatile aroma complexity grape variety", 6),
    ("3.1.1 Methoxypyrazines",
     "3-isobutyl-2-methoxypyrazine IBMP wine green pepper concentration sensory threshold", 8),
    ("3.1.2 Norisoprenoids",
     "beta-damascenone beta-ionone TDN wine norisoprenoid carotenoid degradation", 8),
    ("3.1.3 Fermentation esters",
     "ethyl ester fermentation wine ethyl hexanoate ethyl octanoate yeast", 6),
    ("3.1.4 Monoterpenes",
     "linalool alpha-terpineol monoterpene wine Cabernet glycoside", 6),
    ("3.1.5 Varietal thiols",
     "3-mercaptohexanol 4-mercapto-4-methylpentan-2-one varietal thiol wine", 6),
    ("3.1.6 Brettanomyces ethylphenols",
     "Brettanomyces 4-ethylphenol 4-ethylguaiacol wine fault medicinal", 6),
    ("3.1.7 Smoke taint volatile phenols",
     "smoke taint guaiacol 4-methylguaiacol syringol wine forest fire", 8),
    ("3.1.7b Cork taint TCA",
     "trichloroanisole TCA cork taint wine ng/L detection", 6),
    ("3.2 SPME theory three-phase equilibrium",
     "SPME extraction theory three-phase equilibrium partition coefficient Pawliszyn", 6),
    ("3.3 SPME Arrow geometry capacity",
     "SPME Arrow geometry sorbent volume capacity tube-in-needle DVB CWR PDMS", 8),
    ("3.4 Triple-quadrupole MRM theory",
     "triple quadrupole mass spectrometry MRM scheduled cycle time collision energy", 6),
    ("4.6 SPME extraction parameters",
     "headspace SPME wine extraction temperature time NaCl ethanol optimization", 8),
    ("5.2 DoE optimization Box-Behnken",
     "Box-Behnken response surface design wine SPME headspace optimization", 6),
    ("5.3 Extraction isotherms",
     "SPME extraction isotherm kinetics equilibrium time compound class wine", 6),
    ("6 MRM optimization wine analytes",
     "GC-MS/MS MRM wine volatile collision energy optimization scheduled", 6),
    ("7.1 Linearity validation",
     "matrix-matched calibration linearity R-squared wine analyte SPME", 6),
    ("7.2 LOD LOQ wine analytical",
     "limit of detection quantification LOD LOQ wine GC-MS sensitivity", 6),
    ("7.4 Precision intermediate Horwitz",
     "intermediate precision RSD HorRat repeatability wine validation", 6),
    ("7.5 Matrix effect slope ratio",
     "matrix effect slope ratio wine GC-MS SPME ion suppression", 6),
    ("8 Application multivariate Cab Sauv",
     "Cabernet Sauvignon appellation differentiation PCA OPLS-DA multivariate", 6),
    ("8 Oak compounds Cabernet differentiation",
     "oak vanillin whisky lactone eugenol cis trans wine Cabernet", 6),
]


def fmt_row(r: dict) -> str:
    sim = r.get("boosted", r.get("similarity", 0))
    title = (r.get("title") or "(no title)").strip().replace("\n", " ")
    if len(title) > 180:
        title = title[:177] + "..."
    year = r.get("year") or "?"
    pid = (r.get("paper_id") or "")[:12]
    distill = r.get("distill") or {}
    core = (distill.get("core_idea") or "").strip().replace("\n", " ")
    if len(core) > 250:
        core = core[:247] + "..."
    sig = distill.get("signal_score")
    sig_str = f" [sig {sig}/10]" if sig else ""
    line = f"- **[{sim:.3f}]** {year} · {title} · `{pid}`{sig_str}"
    if core:
        line += f"\n  - tldr: {core}"
    return line


def main() -> int:
    out_path = LIFE_ROOT / "projects" / "wine-spme-dissertation" / "_research" / "rag_synthesis.md"
    lines: list[str] = [
        "# RAG Synthesis — Pgvector Hits per Dissertation Section",
        "",
        "Each block below is a `research search ...` query over the `wine-spme-gcms` "
        "topic (250 papers, nomic-embed-text). Similarity scores include "
        "mark/velocity/distill boosts. Use these as the citation backbone "
        "for the corresponding dissertation section.",
        "",
        "---",
        "",
    ]
    summary: list[str] = []

    for heading, query, k in QUERIES:
        print(f"querying: {heading[:60]}...", flush=True)
        try:
            rows = search(query, k=k, topic="wine-spme-gcms")
        except Exception as e:
            lines.append(f"## §{heading}")
            lines.append("")
            lines.append(f"**ERROR**: `{e!r}`")
            lines.append("")
            continue
        lines.append(f"## §{heading}")
        lines.append(f"*Query: `{query}`*")
        lines.append("")
        if not rows:
            lines.append("(no hits)")
        else:
            for r in rows:
                lines.append(fmt_row(r))
            summary.append(f"- §{heading}: {len(rows)} hits, top score {max(r.get('boosted', r.get('similarity', 0)) for r in rows):.3f}")
        lines.append("")

    # Final summary block
    lines.insert(7, "## Summary")
    lines.insert(8, "")
    for s in summary:
        lines.insert(9, s)
    lines.insert(9 + len(summary), "")
    lines.insert(10 + len(summary), "---")
    lines.insert(11 + len(summary), "")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out_path} ({len(lines)} lines, {len(QUERIES)} sections)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
