# Layout benchmark

How broadly chemcard's 2D layouts work, measured on molecules a chemistry
explanation actually draws, rather than on a hand-picked list of hard cases.

| file | what |
|---|---|
| `molecules.csv` | 231 molecules in 7 categories (drugs, opioids/morphinans, steroids/terpenoids, alkaloids, cages/bridged, macrocycles, teaching examples), structures and CIDs from PubChem |
| `REPORT.md` | latest output of `scripts/layout_report.py` |

```bash
uv run python scripts/layout_report.py --md benchmarks/REPORT.md          # regenerate the report
uv run python scripts/layout_report.py --sheet /tmp/flagged.html          # look at every flagged card
uv run python scripts/fetch_benchmark_molecules.py                        # rebuild the CSV (needs internet)
```

Each molecule goes through the full card pipeline; the final layout is scored with
`render.layout_quality`. "Clean" means no unmarked bond crossing, no clashing atoms,
no atom on a bond, no substituent inside a ring. For each flagged molecule the report
also scores every candidate the layout code already generates: if one was clean, the
selection made a mistake; if none was, a new generator or template is needed.

## Findings (2026-11-19, v0.4 + bridge layouts)

85 % of the set is clean (drugs 99 %, teaching examples 98 %). Of the 35 flagged
molecules, **none is a selection mistake**: no existing candidate was clean. Looked at
one by one (`--sheet`), they fall into four groups:

| group | molecules | verdict |
|---|---|---|
| drawn badly | buprenorphine (improved, see below), huperzine A, picrotoxinin, ginkgolide B, vinblastine | needs new candidates |
| template, crowded label | naloxone, naltrexone, oxycodone (14-OH sits on the bridge) | inherent to the textbook view, could be tuned |
| metric false alarm | adamantanes, cubane, twistane, triptycene (a 3D cage always crosses in 2D); morphine-type templates (the textbook view has one crossing); rifampicin, paclitaxel, clarithromycin, cyclodextrin (substituents pointing into a macrocycle are normal) | drawings are fine |
| out of scope | buckminsterfullerene | no good 2D drawing exists |

Buprenorphine was the one rule problem: the epoxymorphinan template "won outright",
but buprenorphine's extra 6,14-bridge is not part of the template and got drawn
across the rings (6 crossings vs 2 for CoordGen). Fixed the same day: a template now
wins outright only if it adds no crossing beyond its own textbook bridge; otherwise
it competes on score. Buprenorphine now gets CoordGen (2 crossings, substituents
outside the cage); the other nine morphinans keep their templates.

## 2026-11-19 — projected perspectives for any bicyclo[a.b.c]

Clean share unchanged (85 %), no new flags, still no selection mistakes. Twelve
molecules now get a projected perspective view (checked by eye):

- [2.2.2]: bicyclo[2.2.2]octane, DABCO, quinuclidine, barrelene, quinine (its quinuclidine)
- [2.2.1] outside the hand-tuned family: norbornene, norbornadiene, 7-oxanorbornane, epibatidine
- [3.2.1]: tropinone, atropine
- [3.1.1]: beta-pinene

Crowded projections (cocaine, fenchone, 1,8-cineole) lose to the general engines as
intended. Molecules above 40 heavy atoms skip the 3D embedding (amphotericin B took
1.5 s for no usable view); slowest render is back to ~0.3 s.
