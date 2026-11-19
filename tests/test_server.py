"""MCP wiring: tool metadata, UI resource, description budgets."""

import pytest

from chemcard import descriptions
from chemcard.server import CARD_URI, mcp


@pytest.fixture
def anyio_backend():
    return "asyncio"


@pytest.mark.anyio
async def test_tool_is_bound_to_widget():
    tools = {t.name: t for t in await mcp.list_tools()}
    tool = tools["draw_molecule"]
    assert tool.meta["ui"]["resourceUri"] == CARD_URI
    a = tool.annotations
    assert a.read_only_hint is True
    assert a.destructive_hint is False
    assert a.open_world_hint is False
    assert set(tool.input_schema["required"]) == {"smiles"}


@pytest.mark.anyio
async def test_svg_fallback_tool_is_app_only():
    tools = {t.name: t for t in await mcp.list_tools()}
    assert tools["render_structure_svg"].meta["ui"]["visibility"] == ["app"]


def test_svg_fallback_tool_returns_svg():
    from chemcard.server import render_structure_svg

    r = render_structure_svg("CCO")
    assert "<svg" in r.structured_content["svg"]


@pytest.mark.anyio
async def test_widget_resource_mime_type():
    contents = list(await mcp.read_resource(CARD_URI))
    assert contents[0].mime_type == "text/html;profile=mcp-app"
    assert "<!doctype html>" in contents[0].content.lower()


def test_result_keeps_svg_out_of_model_view():
    from chemcard.server import show_structure

    r = show_structure("Cn1cnc2c1c(=O)n(C)c(=O)n2C", label="咖啡因")
    assert not r.is_error
    assert r.structured_content == {
        "kind": "structure",
        "canonical_smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
        "formula": "C8H10N4O2",
        "label": "咖啡因",
    }
    assert "<svg" in r.meta["svg"]
    assert "svg" not in r.content[0].text.lower()


def test_invalid_smiles_is_tool_error():
    from chemcard.server import show_structure

    r = show_structure("C1CC")
    assert r.is_error
    assert "unclosed ring" in r.content[0].text


def test_instructions_fit_chatgpt_budget():
    # ChatGPT relies most on the first ~512 characters of server instructions.
    assert len(descriptions.SERVER_INSTRUCTIONS) <= 512


def test_descriptions_ask_for_highlights():
    # Claude test 2026-11-19: highlights existed but were never used until the text asked for them
    assert "highlights" in descriptions.SERVER_INSTRUCTIONS
    assert "highlight" in descriptions.SHOW_STRUCTURE and "highlights" in descriptions.SHOW_REACTION_STEP


@pytest.mark.anyio
async def test_schemas_have_no_refs_for_claude():
    # Claude drops $ref: nested fields (Species.highlights, MechanismStep.*) became `items: {}`
    import json

    tools = {t.name: t for t in await mcp.list_tools()}
    for t in tools.values():
        dumped = json.dumps(t.input_schema)
        assert "$ref" not in dumped and "$defs" not in dumped, t.name
    step_props = tools["draw_reaction_step"].input_schema["properties"]
    assert {"smiles", "label"} <= step_props["reactants"]["items"]["properties"].keys()
    assert {"side", "index", "smarts", "atoms", "color"} <= step_props["highlights"]["items"]["properties"].keys()
    step = tools["check_mechanism"].input_schema["properties"]["steps"]["items"]["properties"]
    assert {"reactants", "products", "conditions", "radical"} <= step.keys()


@pytest.mark.anyio
async def test_no_schema_deeper_than_show_structure():
    # show_structure's schema is confirmed to reach the model intact in Claude; keep others no deeper
    def depth(n, d=0):
        if isinstance(n, dict):
            return max([depth(v, d + 1) for v in n.values()] + [d])
        if isinstance(n, list):
            return max([depth(v, d + 1) for v in n] + [d])
        return d

    tools = {t.name: t for t in await mcp.list_tools()}
    limit = depth(tools["draw_molecule"].input_schema)
    for t in tools.values():
        assert depth(t.input_schema) <= limit, t.name
