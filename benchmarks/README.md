# Layout benchmark

How broadly chemcard's 2D layouts work, measured on molecules a chemistry
explanation actually draws, rather than on a hand-picked list of hard cases.

| file | what |
|---|---|
| `molecules.csv` | 231 samples in 7 categories (drugs, opioids/morphinans, steroids/terpenoids, alkaloids, cages/bridged, macrocycles, teaching examples) = **230 unique compounds** (eucalyptol and 1,8-cineole are both CID 2758); structures and CIDs from PubChem |
| `REPORT.md` | latest output of `scripts/layout_report.py` |

```bash
uv run python scripts/layout_report.py --md benchmarks/REPORT.md          # regenerate the report
uv run python scripts/layout_report.py --sheet /tmp/flagged.html          # look at every flagged card
uv run python scripts/fetch_benchmark_molecules.py                        # rebuild the CSV (needs internet)
```

`fetch_benchmark_molecules.py` retries timeouts and 429/5xx; if any name still fails
it lists the failures, exits non-zero and leaves `molecules.csv` untouched. Only a
complete result replaces the file, atomically.

## What the report measures

Each sample goes through the full card pipeline. Two kinds of problems are counted:

- **layout problems**, from the chosen coordinates (`render.layout_quality`): an
  unmarked bond crossing, clashing atoms, an atom on a bond, a substituent inside a ring;
- **render overlaps**, from the final SVG (`render.svg_label_overlaps`): two labels
  overlapping, or a bond running through a label it is not attached to. These catch
  atoms RDKit adds only while drawing — in artemisinin an explicit H on a
  stereocentre sat 3.7 px from the peroxide O while the coordinates looked clean.

"Clean" means neither. Category rows count samples; the overall rate counts each
compound (CID) once.

For each flagged sample the report lays out the other candidates and shows the
trade-off against the chosen one, e.g. `coordgen: crossing −4, clash +1`. Its summary
claims only what those counts support:

- **missed fully clean candidates** — a candidate with all four layout counts at zero
  existed but was not chosen;
- **Pareto-better candidates** — another candidate is no worse on any layout count and
  better on one;
- **trade-offs** — one candidate is better than the chosen drawing on some count and
  worse on another. Candidates are judged one at a time, so a molecule can have both a
  Pareto-better candidate and a trade-off.

None of these says a choice was right or wrong. "Zero missed clean candidates" only
means no flawless alternative was overlooked: before the template fix, buprenorphine
was drawn with the template (6 crossings) instead of CoordGen (2 crossings, 1 more
clash), and that count was still zero. Render overlaps are not compared across
candidates (that would need every candidate drawn).

## Findings (2026-11-19)

83 % of the 230 unique compounds are clean (drugs 99 %, teaching examples 98 %).
16 samples have render overlaps; four of them had looked clean before this check
(artemisinin, forskolin, epibatidine, tacrolimus — all checked by eye: a drawing-time
H on an O, or a bond through an OH/NH label).

Missed fully clean candidates: 0. Pareto-better candidates: 4 — dextromethorphan,
naloxone, naltrexone and oxycodone, where RDKit's layout has fewer layout problems
than the morphinan template. The template is kept on purpose: the textbook view is
the convention, and RDKit's alternative stretches the ether bridge across the drawing
(compared by eye). This is exactly the kind of judgement the counts cannot make.

Looked at one by one (`--sheet`), the flagged samples fall into these groups:

| group | molecules | verdict |
|---|---|---|
| drawn badly | buprenorphine (improved, see below), huperzine A, picrotoxinin, ginkgolide B, vinblastine | needs new candidates |
| render overlaps | artemisinin, forskolin, epibatidine, tacrolimus, and on top of other problems in several morphinans, picrotoxinin, ginkgolide B, cyclodextrin and three macrocycles | label placement; not addressed yet |
| template, crowded label | naloxone, naltrexone, oxycodone (14-OH sits on the bridge) | inherent to the textbook view, could be tuned |
| metric false alarm | adamantanes, cubane, twistane, triptycene (a 3D cage always crosses in 2D); morphine-type templates (the textbook view has one crossing); rifampicin, paclitaxel, clarithromycin, cyclodextrin (substituents pointing into a macrocycle are normal) | drawings are fine |
| out of scope | buckminsterfullerene | no good 2D drawing exists |

Buprenorphine was a rule problem: the epoxymorphinan template "won outright", but
buprenorphine's extra 6,14-bridge is not part of the template and got drawn across
the rings (6 crossings vs 2 for CoordGen). Fixed: a template now wins outright only if
it adds no crossing beyond its own textbook bridge; otherwise it competes on score.
Buprenorphine now gets CoordGen; the other nine morphinans keep their templates.
