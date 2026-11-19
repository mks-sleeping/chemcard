"""Text the model reads to decide *when* and *how* to call the tools.

Kept apart from the code because it is tuned repeatedly against evals/.
Anything the model needs in order to call a tool correctly must live here
(tool or parameter descriptions), never only in an optional Skill.
"""

# Claude loads connector tools on demand: before a tool is loaded the model only sees short
# text (server name/title/description, tool names). 2026-11-19, fresh chat: "阿司匹林哪部分会水解？"
# was answered with no tool at all. So the server text says plainly WHEN to draw.
SERVER_DESCRIPTION = (
    "Draws clean 2D chemical structures, reaction steps and mechanisms inline in the answer. "
    "Use it whenever you explain a specific molecule or reaction."
)

# ChatGPT weighs roughly the first 512 characters; keep the essentials first.
SERVER_INSTRUCTIONS = (
    "Whenever an answer is about a specific molecule or reaction (its structure, which part reacts, "
    "a mechanism), draw it inline like a textbook figure: draw_molecule for a molecule, "
    "draw_reaction_step for each step (check_mechanism first if multi-step). Interleave: figure, "
    "then its explanation. Mark the atoms/bonds you discuss with `highlights`. If the user selected "
    "atoms (mapped SMILES in their message), redraw from it and highlight them. Skip trivial species."
)

SHOW_STRUCTURE = (
    "Show a clean 2D structure of ONE molecule inline in the explanation. "
    "Call it when the explanation refers to a specific molecule's structure "
    "(functional groups, protecting groups, stereochemistry, why it reacts as it does). "
    "Put each figure next to the paragraph about that molecule: when comparing two molecules, "
    "draw the first, explain it, then draw the second and explain it; do not draw both first. "
    "When the text is about part of the molecule (a group, an atom, a bond), highlight that part "
    "with `highlights` so the reader sees what the sentence refers to. If the user selected atoms "
    "in a card, redraw from the atom-mapped SMILES they sent and highlight those atom numbers. "
    "For a reaction or mechanism step use draw_reaction_step instead. "
    "Do not call it for trivial species (H2O, HCl) or purely conceptual answers. "
    "If the tool reports invalid SMILES, fix the SMILES and call again."
)

SHOW_REACTION_STEP = (
    "Show ONE reaction or mechanism step as a scheme 'A + B → C' (conditions over the arrow). "
    "When explaining a mechanism, go step by step: call this tool for step 1, then explain step 1 "
    "directly below it in 1–4 sentences, then call it for step 2, and so on. Never draw all steps "
    "first and explain them afterwards. A step may have 1–4 species per side. Small byproducts "
    "(H2O, H+) can go in the products or in `conditions`. Mark the reacting centres with `highlights` "
    "(e.g. the attacking atom in blue, the electrophilic carbon or the bond that breaks in red). Write charged atoms with their H explicitly "
    "([NH+], [OH2+]); a bracket atom without H becomes a radical and the step is refused. The result "
    "reports imbalance and unusual species (methyl/primary cations, strong base in acid); if one "
    "reveals a mistake, fix the step and call again."
)

CHECK_MECHANISM = (
    "Validate a whole multi-step mechanism before drawing it (no figure is shown). Checks every "
    "species (valid SMILES, unexpected radicals from a bracket atom missing its H, unstabilised "
    "methyl/primary cations, strong bases under acidic conditions and vice versa, doubly charged "
    "species), each step's atom and charge balance, and that every step starts from an intermediate "
    "produced by the previous step. Call it once with all steps, fix every [error], consider each "
    "[warning], then draw the steps one at a time with draw_reaction_step."
)

# Appended to every successful drawing result. The model reads it right after each
# figure, which steers the "figure → its explanation → next figure" rhythm far more
# reliably than the tool description alone (Claude test 2026-11-19: steps 1–3 were
# drawn back to back without text in between).
NEXT_AFTER_STRUCTURE = (
    "Next: write the explanation of this molecule now, directly below the figure, "
    "before drawing anything else."
)
# Appended when a step was drawn without any highlight. Like NEXT_AFTER_*, feedback right
# after the call steers better than the description (2026-11-19: three tests in a row drew
# steps without highlights although the description asked for them).
NO_HIGHLIGHT_HINT = (
    "Tip: this step has no highlights. For the next steps, pass `highlights` to colour the "
    "reacting atoms/bonds (e.g. {\"side\": \"reactants\", \"index\": 1, \"smarts\": \"C=O\", \"color\": \"red\"})."
)
NEXT_AFTER_STEP = (
    "Next: explain this step now, in 1–4 sentences directly below the figure. "
    "Only after that, draw the next step."
)

SMILES_PARAM = (
    "SMILES of the molecule. Use '*' for a generic substituent R. "
    "Write stereochemistry (@, @@, /, \\) only when it matters for the explanation."
)

LABEL_PARAM = "Short caption shown under the structure, in the user's language (e.g. 'R–O–TBS 硅醚')."

SPECIES_LABEL_PARAM = "Optional short name under this structure, in the user's language (e.g. '亚胺鎓离子 ①')."

CONDITIONS_PARAM = "Reagents/conditions written over the arrow, e.g. 'H⁺, −H₂O' or 'NaOH, H₂O, Δ'. Use Unicode sub/superscripts."

CAPTION_PARAM = "Step title shown under the scheme, in the user's language (e.g. '第 1 步：形成亚胺鎓离子')."

REVERSIBLE_PARAM = "True to draw an equilibrium arrow (⇌) instead of →."

HIGHLIGHTS_PARAM = (
    "Optional coloured halos that point at what the explanation talks about, e.g. "
    "[{\"smarts\": \"C=O\", \"color\": \"red\"}] for the carbonyl, or "
    "[{\"atoms\": [1], \"color\": \"blue\"}] for the atom written [NH2:1] in the SMILES. "
    "Use 1–2 colours; highlight only what the next sentences discuss. Atom-map numbers are not drawn."
)

HIGHLIGHT_SMARTS_PARAM = "SMARTS of the group to highlight (all matches are highlighted), e.g. 'C=O', '[OX2H]', 'c1ccccc1'."

HIGHLIGHT_ATOMS_PARAM = (
    "Atom-map numbers to highlight; write them in the SMILES first, e.g. 'C[C@H]([NH2:1])C(=O)O'. "
    "When the user selected atoms in a card, reuse the numbers from the mapped SMILES they sent."
)

STEP_HIGHLIGHTS_PARAM = (
    "Optional coloured halos on the species of this step, each naming its species by side and "
    "1-based index, e.g. [{\"side\": \"reactants\", \"index\": 1, \"smarts\": \"C(=O)O\", \"color\": \"red\"}, "
    "{\"side\": \"reactants\", \"index\": 2, \"atoms\": [1], \"color\": \"blue\"}] (the second needs "
    "[OH2:1] in that SMILES). Mark what reacts: the attacking atom, the electrophilic atom, the bond that breaks or forms."
)

STEP_HIGHLIGHT_SIDE_PARAM = "Which side the highlighted species is on."

STEP_HIGHLIGHT_INDEX_PARAM = "1-based position of the species on that side (1 = first reactant/product)."

RADICAL_PARAM = "True only for radical chemistry (homolysis, radical chain steps); allows species with unpaired electrons."
