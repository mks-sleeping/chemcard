"""MCP adapter: registers tools and the UI resource, assembles results.

No chemistry here — that lives in `chemcard.chem`.

Result layout (who sees what):
- content            model + widget; the only thing text-only hosts show
- structuredContent  model + widget; keep it small (formula, canonical SMILES, label)
- _meta              widget only; carries the SVGs so they never enter the model context
                     (Claude does not forward _meta to the card; the card then fetches
                     each SVG through the app-only tool render_structure_svg)

One widget (card.html) renders every tool; structuredContent.kind selects the layout.
"""

from __future__ import annotations

import argparse
import os
from importlib.resources import files
from typing import Annotated, Literal

from mcp.server import MCPServer
from mcp.server.apps import Apps
from mcp.server.transport_security import TransportSecuritySettings
from mcp.types import CallToolResult, TextContent, ToolAnnotations
from pydantic import BaseModel, Field

from chemcard import __version__, descriptions
from chemcard.chem import mechanism
from chemcard.chem.highlight import HighlightSpec, InvalidHighlight, describe_selection, resolve
from chemcard.chem.lint import lint_species
from chemcard.chem.render import InvalidSmiles, Structure, render_smiles
from chemcard.chem.scheme import ArrowSpec, InvalidArrow, compose, resolve_arrows

# The URI is a cache key for hosts: bump the version on breaking widget changes.
CARD_URI = "ui://chemcard/card-v6.html"
MAX_SPECIES_PER_SIDE = 4

READ_ONLY = ToolAnnotations(
    read_only_hint=True,
    destructive_hint=False,
    idempotent_hint=True,
    open_world_hint=False,
)

apps = Apps()
apps.add_html_resource(
    CARD_URI,
    (files("chemcard") / "widgets" / "card.html").read_text(encoding="utf-8"),
    title="Chemical structure",
    prefers_border=False,  # look like an inline figure, not a boxed app
)


class Highlight(BaseModel):
    smarts: str | None = Field(default=None, description=descriptions.HIGHLIGHT_SMARTS_PARAM)
    atoms: list[int] | None = Field(default=None, description=descriptions.HIGHLIGHT_ATOMS_PARAM)
    color: Literal["red", "blue", "green", "orange", "purple"] = "red"


def _error(text: str) -> CallToolResult:
    return CallToolResult(content=[TextContent(type="text", text=text)], is_error=True)


def _draw(smiles: str, highlights: list[Highlight] | None) -> tuple[Structure, dict | None]:
    """Resolve highlights, draw from canonical SMILES. Raises InvalidSmiles / InvalidHighlight."""
    specs = [HighlightSpec(color=h.color, smarts=h.smarts, atoms=tuple(h.atoms or ())) for h in highlights or []]
    r = resolve(smiles, specs)
    s = render_smiles(r.canonical_smiles, r.atom_colors, r.bond_colors)
    marks = None
    if r.atom_colors or r.bond_colors:  # small; lets the card redraw identically via the fallback tool
        marks = {"atoms": {str(i): c for i, c in r.atom_colors.items()},
                 "bonds": {str(i): c for i, c in r.bond_colors.items()}}
    return s, marks


def _figure(s: Structure) -> dict:
    """Widget-only payload for one drawing."""
    return {"svg": s.svg, "width": s.width, "height": s.height, "atoms": list(s.atoms)}


@apps.tool(
    resource_uri=CARD_URI,
    name="draw_molecule",
    title="Draw molecule",
    description=descriptions.SHOW_STRUCTURE,
    annotations=READ_ONLY,
)
def show_structure(
    smiles: Annotated[str, Field(description=descriptions.SMILES_PARAM)],
    label: Annotated[str | None, Field(description=descriptions.LABEL_PARAM)] = None,
    highlights: Annotated[list[Highlight] | None, Field(description=descriptions.HIGHLIGHTS_PARAM)] = None,
) -> CallToolResult:
    try:
        s, marks = _draw(smiles, highlights)
    except InvalidSmiles as e:
        return _error(f"Invalid SMILES: {e}")
    except InvalidHighlight as e:
        return _error(f"Not drawn, highlight problem: {e}")
    shown = label or s.formula
    text = f"Displayed {shown} ({s.formula})."
    issues = lint_species(s.canonical_smiles)
    if issues:  # advice only: a single molecule may legitimately be a radical or a cation
        text += " Check: " + " ".join(i.message for i in issues)
    text += "\n" + descriptions.NEXT_AFTER_STRUCTURE
    data = {"kind": "structure", "canonical_smiles": s.canonical_smiles, "formula": s.formula, "label": label}
    if marks:
        data["highlight"] = marks
    return CallToolResult(
        content=[TextContent(type="text", text=text)],
        structured_content=data,
        meta=_figure(s),
    )


class Species(BaseModel):
    smiles: str = Field(description=descriptions.SMILES_PARAM)
    label: str | None = Field(default=None, description=descriptions.SPECIES_LABEL_PARAM)


class StepHighlight(Highlight):
    """A highlight on one species of a step. Kept top-level (not inside Species) so the
    schema stays as shallow as show_structure's: Claude did not show the model a
    Species-level `highlights` field (2026-11-19)."""

    side: Literal["reactants", "products"] = Field(description=descriptions.STEP_HIGHLIGHT_SIDE_PARAM)
    index: int = Field(ge=1, le=MAX_SPECIES_PER_SIDE, description=descriptions.STEP_HIGHLIGHT_INDEX_PARAM)


class Arrow(BaseModel):
    start: str = Field(description=descriptions.ARROW_START_PARAM)
    end: str = Field(description=descriptions.ARROW_END_PARAM)
    style: Literal["full", "half"] = Field(default="full", description=descriptions.ARROW_STYLE_PARAM)


def _marks_to_colors(marks: dict | None) -> tuple[dict[int, str], dict[int, str]]:
    m = marks or {}
    return ({int(i): c for i, c in (m.get("atoms") or {}).items()},
            {int(i): c for i, c in (m.get("bonds") or {}).items()})


def _scheme_figure(sc) -> dict:
    return {"svg": sc.svg, "width": sc.width, "height": sc.height, "atoms": sc.atoms, "species": sc.species}


@apps.tool(
    resource_uri=CARD_URI,
    name="draw_reaction_step",
    title="Draw reaction step",
    description=descriptions.SHOW_REACTION_STEP,
    annotations=READ_ONLY,
)
def show_reaction_step(
    reactants: Annotated[list[Species], Field(min_length=1, max_length=MAX_SPECIES_PER_SIDE)],
    products: Annotated[list[Species], Field(min_length=1, max_length=MAX_SPECIES_PER_SIDE)],
    conditions: Annotated[str | None, Field(description=descriptions.CONDITIONS_PARAM)] = None,
    caption: Annotated[str | None, Field(description=descriptions.CAPTION_PARAM)] = None,
    reversible: Annotated[bool, Field(description=descriptions.REVERSIBLE_PARAM)] = False,
    radical: Annotated[bool, Field(description=descriptions.RADICAL_PARAM)] = False,
    highlights: Annotated[list[StepHighlight] | None, Field(description=descriptions.STEP_HIGHLIGHTS_PARAM)] = None,
    arrows: Annotated[list[Arrow] | None, Field(description=descriptions.ARROWS_PARAM)] = None,
) -> CallToolResult:
    sides = {"reactants": reactants, "products": products}
    for h in highlights or []:
        if h.index > len(sides[h.side]):
            return _error(f"highlight points at {h.side}[{h.index}], but there are only {len(sides[h.side])} {h.side}")

    step = mechanism.Step(
        reactants=[sp.smiles for sp in reactants],
        products=[sp.smiles for sp in products],
        conditions=conditions,
        radical=radical,
    )
    try:
        issues = mechanism.lint_step(step)
    except InvalidSmiles as e:
        return _error(f"Invalid SMILES in {e}")
    errors = [i for i in issues if i.level == "error"]
    if errors:  # almost certainly a typo: refuse to draw so the model fixes it
        return _error("Not drawn. Fix these and call again:\n" + "\n".join(i.text() for i in errors))
    warnings = [i.text() for i in issues]

    drawn: dict[str, list[tuple[Structure, dict | None]]] = {}
    for side, species in (("reactants", reactants), ("products", products)):
        drawn[side] = []
        for i, sp in enumerate(species, start=1):
            try:
                mine = [h for h in highlights or [] if h.side == side and h.index == i]
                drawn[side].append(_draw(sp.smiles, mine))
            except InvalidSmiles as e:
                return _error(f"Invalid SMILES in {side}[{i}]: {e}")
            except InvalidHighlight as e:
                return _error(f"Not drawn, highlight problem in {side}[{i}]: {e}")

    resolved_arrows: list[dict] = []
    reactant_scheme = None
    if arrows:
        try:
            canonicals, resolved_arrows = resolve_arrows(
                [sp.smiles for sp in reactants],
                [ArrowSpec(src=a.start, dst=a.end, style=a.style) for a in arrows],
            )
            reactant_scheme = compose(canonicals, resolved_arrows, [_marks_to_colors(m) for _, m in drawn["reactants"]])
        except InvalidArrow as e:
            return _error(f"Not drawn, arrow problem: {e}")

    arrow = " ⇌ " if reversible else " → "
    scheme = arrow.join(
        " + ".join(sp.label or st.formula for sp, (st, _) in zip(species, drawn[side]))
        for side, species in (("reactants", reactants), ("products", products))
    )
    text = f"Displayed step: {scheme}."
    if arrows:
        text += " Curved arrows: " + "; ".join(f"{a.start} → {a.end}" + (" (half)" if a.style == "half" else "") for a in arrows) + "."
    if warnings:
        text += "\nCheck (drawn anyway; fix and redraw if it is a mistake):\n" + "\n".join(warnings)
    if not highlights:
        text += "\n" + descriptions.NO_HIGHLIGHT_HINT
    text += "\n" + descriptions.NEXT_AFTER_STEP

    def summary(species: list[Species], items: list[tuple[Structure, dict | None]]) -> list[dict]:
        out = []
        for sp, (st, marks) in zip(species, items):
            d = {"label": sp.label, "canonical_smiles": st.canonical_smiles, "formula": st.formula}
            if marks:
                d["highlight"] = marks
            out.append(d)
        return out

    return CallToolResult(
        content=[TextContent(type="text", text=text)],
        structured_content={
            "kind": "reaction",
            "caption": caption,
            "conditions": conditions,
            "reversible": reversible,
            "reactants": summary(reactants, drawn["reactants"]),
            "products": summary(products, drawn["products"]),
            "warnings": warnings,
            **({"arrows": resolved_arrows} if resolved_arrows else {}),
        },
        meta={
            **{side: [_figure(st) for st, _ in drawn[side]] for side in ("reactants", "products")},
            **({"reactant_scheme": _scheme_figure(reactant_scheme)} if reactant_scheme else {}),
        },
    )


@apps.tool(
    resource_uri=CARD_URI,
    visibility=["app"],
    name="render_reactant_scheme",
    title="Render reactants with arrows",
    description="Card-internal: redraw a step's reactants with curved arrows when _meta was not delivered.",
    annotations=READ_ONLY,
)
def render_reactant_scheme(reactants: list[dict], arrows: list[dict]) -> CallToolResult:
    """`reactants`: [{"smiles": canonical, "highlight": {...} | None}]; `arrows`: resolved, from structuredContent."""
    try:
        sc = compose(
            [r["smiles"] for r in reactants],
            arrows,
            [_marks_to_colors(r.get("highlight")) for r in reactants],
        )
    except (InvalidSmiles, InvalidArrow, KeyError, IndexError, TypeError, ValueError) as e:
        return _error(str(e))
    return CallToolResult(content=[TextContent(type="text", text="ok")], structured_content=_scheme_figure(sc))


@apps.tool(
    resource_uri=CARD_URI,
    visibility=["app"],  # hidden from the model; only the card calls it
    name="render_structure_svg",
    title="Render structure SVG",
    description="Card-internal: fetch a drawing when the host did not deliver it in the tool result's _meta.",
    annotations=READ_ONLY,
)
def render_structure_svg(smiles: str, highlight: dict | None = None) -> CallToolResult:
    """`smiles` is canonical; `highlight` is the {"atoms": {...}, "bonds": {...}} from structuredContent."""
    h = highlight or {}
    try:
        s = render_smiles(
            smiles,
            {int(i): c for i, c in (h.get("atoms") or {}).items()},
            {int(i): c for i, c in (h.get("bonds") or {}).items()},
        )
    except (InvalidSmiles, KeyError, ValueError) as e:
        return _error(str(e))
    return CallToolResult(content=[TextContent(type="text", text="ok")], structured_content=_figure(s))


@apps.tool(
    resource_uri=CARD_URI,
    visibility=["app"],
    name="describe_selection",
    title="Describe selected atoms",
    description="Card-internal: turn atoms the user clicked into fragment and atom-mapped SMILES.",
    annotations=READ_ONLY,
)
def describe_selection_tool(smiles: str, atoms: list[int]) -> CallToolResult:
    try:
        d = describe_selection(smiles, atoms)
    except InvalidSmiles as e:
        return _error(str(e))
    return CallToolResult(content=[TextContent(type="text", text=d["fragment_smiles"])], structured_content=d)


mcp = MCPServer(
    name="chemcard",
    title="Chemistry structure drawings",
    description=descriptions.SERVER_DESCRIPTION,
    version=__version__,
    instructions=descriptions.SERVER_INSTRUCTIONS,
    extensions=[apps],
)


class MechanismStep(BaseModel):
    reactants: list[str] = Field(min_length=1, max_length=MAX_SPECIES_PER_SIDE, description="SMILES of the species entering this step.")
    products: list[str] = Field(min_length=1, max_length=MAX_SPECIES_PER_SIDE, description="SMILES of the species leaving this step.")
    conditions: str | None = Field(default=None, description=descriptions.CONDITIONS_PARAM)
    radical: bool = Field(default=False, description=descriptions.RADICAL_PARAM)


@mcp.tool(
    name="check_mechanism",
    title="Check mechanism",
    description=descriptions.CHECK_MECHANISM,
    annotations=READ_ONLY,
)
def check_mechanism(
    steps: Annotated[list[MechanismStep], Field(min_length=1, max_length=12)],
) -> CallToolResult:
    report = mechanism.check_mechanism([mechanism.Step(**st.model_dump()) for st in steps])
    return CallToolResult(
        content=[TextContent(type="text", text=report.text())],
        structured_content=report.as_dict(),
    )


def flatten_schema(schema: dict) -> dict:
    """Inline `$ref`/`$defs` and turn `anyOf: [X, null]` into X.

    Claude's client drops `$ref`: nested models (Species, Highlight,
    MechanismStep) reached the model as `items: {}`, so it never saw
    `highlights` (found 2026-11-19). Pydantic still validates the real
    arguments; this only changes what the model is shown.
    """
    defs = schema.get("$defs", {})

    def walk(node):
        if isinstance(node, list):
            return [walk(n) for n in node]
        if not isinstance(node, dict):
            return node
        if "$ref" in node:
            target = walk(defs[node["$ref"].split("/")[-1]])
            extra = {k: walk(v) for k, v in node.items() if k != "$ref"}
            return {**target, **extra}
        out = {k: walk(v) for k, v in node.items() if k not in ("$defs", "title")}
        options = out.get("anyOf")
        if isinstance(options, list):
            non_null = [o for o in options if o != {"type": "null"}]
            if len(non_null) == 1 and len(non_null) < len(options):
                out.pop("anyOf")
                out = {**non_null[0], **out}
        return out

    return walk(schema)


# Applied once every tool is registered (uses the SDK's tool registry, which has no public setter).
for _tool in mcp._tool_manager.list_tools():
    _tool.parameters = flatten_schema(_tool.parameters)


def _transport_security() -> TransportSecuritySettings | None:
    """Allow tunnel hostnames (ngrok/cloudflared) listed in CHEMCARD_ALLOWED_HOSTS.

    Unset: the SDK's default DNS-rebinding protection accepts localhost only,
    and a tunnel request fails with `421 Invalid Host header`.
    """
    raw = os.environ.get("CHEMCARD_ALLOWED_HOSTS", "")
    hosts = [h.strip() for h in raw.split(",") if h.strip()]
    if not hosts:
        return None
    local = ["127.0.0.1:*", "localhost:*", "[::1]:*"]
    return TransportSecuritySettings(
        enable_dns_rebinding_protection=True,
        allowed_hosts=local + [x for h in hosts for x in (h, f"{h}:*")],
        allowed_origins=[f"https://{h}" for h in hosts],
    )


def main() -> None:
    parser = argparse.ArgumentParser(prog="chemcard")
    parser.add_argument("--http", action="store_true", help="serve Streamable HTTP at /mcp instead of stdio")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()

    if not args.http:
        mcp.run("stdio")
        return

    import uvicorn

    app = mcp.streamable_http_app(
        stateless_http=True,
        json_response=True,
        transport_security=_transport_security(),
        host=args.host,
    )
    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
