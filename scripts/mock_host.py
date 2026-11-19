"""Local mock host: preview the cards without Claude or ChatGPT.

Serves a page that embeds the real widget in sandboxed iframes and plays the
host side of MCP Apps: ui/initialize → initialized → tool-input → tool-result,
answers the card's tools/call by calling the server functions in this process,
and logs ui/message / ui/update-model-context. Text paragraphs sit between the
cards like a chat answer.

    uv run python scripts/mock_host.py [--port 8765]

then open http://127.0.0.1:8765/ (all samples, with message logs) or
http://127.0.0.1:8765/?demo=1 (a clean chat-like page, used for the README image).
This checks our own bridge against our reading of the protocol; it does not
replace testing in a real host.
"""

from __future__ import annotations

import argparse
import html
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import anyio

from chemcard.server import (
    CARD_URI,
    Arrow,
    Highlight,
    Species,
    StepHighlight,
    describe_selection_tool,
    mcp,
    render_reactant_scheme,
    render_structure_svg,
    show_reaction_step,
    show_structure,
)

S, H = Species, Highlight

# (paragraph, tool, arguments, host behaviour)
#   normal      host forwards the full CallToolResult
#   strip_meta  host drops _meta (as Claude does); the card fetches drawings via tools/call
#   silent      host never sends tool-result; the card must show its diagnostics
SAMPLES = [
    ("阿司匹林里有两个羰基：红色是酯（会水解的部分），蓝色是羧酸。点任意原子试试选中：", "show_structure",
     {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "label": "阿司匹林",
      "highlights": [H(smarts="[CX3](=O)[OX2][#6]", color="red"), H(smarts="C(=O)[OX2H1]", color="blue")]}, "strip_meta"),
    ("用原子编号指定单个原子：L-丙氨酸的氨基氮。", "show_structure",
     {"smiles": "C[C@H]([NH2:1])C(=O)O", "label": "L-丙氨酸", "highlights": [H(atoms=[1], color="green")]}, "normal"),
    ("机理分步展示：亲核的氮进攻醛基（橙色）。", "show_reaction_step",
     {"reactants": [S(smiles="O=CCCC=O", label="丁二醛"), S(smiles="C[NH2:1]", label="甲胺")],
      "products": [S(smiles="CN1C(O)CCC1O", label="2,5-二羟基-1-甲基吡咯烷")],
      "conditions": "H₂O", "caption": "第 1 步：两次加成，关环成吡咯烷",
      "highlights": [StepHighlight(side="reactants", index=1, smarts="[CX3H1]=O", color="orange"),
                     StepHighlight(side="reactants", index=2, atoms=[1], color="blue")]}, "strip_meta"),
    ("酸催化下脱去一分子水，得到亲电的亚胺鎓离子：", "show_reaction_step",
     {"reactants": [S(smiles="CN1C(O)CCC1O"), S(smiles="[H+]")],
      "products": [S(smiles="C[N+]1=CCCC1O", label="亚胺鎓离子 ①"), S(smiles="O")],
      "conditions": "H⁺", "caption": "第 2 步：脱水生成亚胺鎓离子", "reversible": True,
      "highlights": [StepHighlight(side="products", index=1, smarts="C=[N+]", color="red")]}, "normal"),
    ("最后加热脱去两分子 CO₂：", "show_reaction_step",
     {"reactants": [S(smiles="CN1C2CCC1C(C(=O)O)C(=O)C2C(=O)O", label="托品酮-2,4-二羧酸")],
      "products": [S(smiles="CN1C2CCC1CC(=O)C2", label="托品酮"), S(smiles="O=C=O"), S(smiles="O=C=O")],
      "conditions": "Δ", "caption": "第 6 步：两次脱羧"}, "normal"),
    ("SN2：OH⁻ 的孤对电子从背面进攻碳，C–Br 键的电子流向 Br（弯箭头）：", "show_reaction_step",
     {"reactants": [S(smiles="[OH-:1]", label="氢氧根"), S(smiles="[CH3:2][Br:3]", label="溴甲烷")],
      "products": [S(smiles="CO", label="甲醇"), S(smiles="[Br-]")],
      "caption": "SN2：一步协同", "arrows": [Arrow(start="lp:1", end="bond:1-2"), Arrow(start="bond:2-3", end="atom:3")]},
     "normal"),
    ("氰根加成到酮羰基（宿主丢弃 _meta，卡片自己取图）：", "show_reaction_step",
     {"reactants": [S(smiles="N#[C-:1]", label="氰根"), S(smiles="C[C:2](C)=[O:3]", label="丙酮")],
      "products": [S(smiles="CC(C)(C#N)[O-]", label="氰醇负离子")],
      "caption": "亲核加成", "highlights": [StepHighlight(side="reactants", index=2, atoms=[2], color="red")],
      "arrows": [Arrow(start="lp:1", end="atom:2"), Arrow(start="bond:2-3", end="atom:3")]}, "strip_meta"),
    ("E2：碱夺取 β-H，C–H 电子形成 C=C，Br⁻ 离去：", "show_reaction_step",
     {"reactants": [S(smiles="CC[O-:1]", label="乙氧负离子"), S(smiles="C[CH:2]([H:5])[CH2:3][Br:4]", label="1-溴丙烷")],
      "products": [S(smiles="CC=C"), S(smiles="CCO"), S(smiles="[Br-]")],
      "caption": "E2", "arrows": [Arrow(start="lp:1", end="atom:5"), Arrow(start="bond:2-5", end="bond:2-3"),
                                  Arrow(start="bond:3-4", end="atom:4")]}, "normal"),
    ("自由基均裂：两个单电子箭头（鱼钩）：", "show_reaction_step",
     {"reactants": [S(smiles="[Br:1][Br:2]")], "products": [S(smiles="[Br]"), S(smiles="[Br]")],
      "conditions": "hν", "radical": True, "caption": "均裂",
      "arrows": [Arrow(start="bond:1-2", end="atom:1", style="half"), Arrow(start="bond:1-2", end="atom:2", style="half")]},
     "normal"),
    ("错误的 SMILES 不显示卡片，模型会收到错误信息后重试：", "show_structure", {"smiles": "C1CC"}, "normal"),
    ("【宿主不发结果】4 秒后卡片应显示调试信息：", "show_structure", {"smiles": "CCO", "label": "乙醇"}, "silent"),
]

# Clean page for the README screenshot: explanation → figure → explanation.
DEMO = [
    ("阿司匹林会水解的是酯基（红色）；同一分子里的羧基（蓝色）本身不会水解。", "show_structure",
     {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "label": "阿司匹林",
      "highlights": [H(smarts="[C;$(C(=O)Oc)](=O)O", color="red"), H(smarts="C(=O)[OX2H1]", color="blue")]}, "normal"),
    ("SN2 一步完成：OH⁻ 的孤对电子从背面进攻碳，C–Br 键的电子流向 Br：", "show_reaction_step",
     {"reactants": [S(smiles="[OH-:1]", label="氢氧根"), S(smiles="[CH3:2][Br:3]", label="溴甲烷")],
      "products": [S(smiles="CO", label="甲醇"), S(smiles="[Br-]")], "caption": "SN2",
      "highlights": [StepHighlight(side="reactants", index=2, atoms=[2], color="red")],
      "arrows": [Arrow(start="lp:1", end="bond:1-2"), Arrow(start="bond:2-3", end="atom:3")]}, "normal"),
    ("桥环骨架也能画清楚，例如托品酮合成的最后一步：", "show_reaction_step",
     {"reactants": [S(smiles="CN1C2CCC1C(C(=O)O)C(=O)C2C(=O)O", label="托品酮-2,4-二羧酸")],
      "products": [S(smiles="CN1C2CCC1CC(=O)C2", label="托品酮"), S(smiles="O=C=O"), S(smiles="O=C=O")],
      "conditions": "Δ", "caption": "两次脱羧"}, "normal"),
]

DRAW_TOOLS = {"show_structure": show_structure, "show_reaction_step": show_reaction_step}
APP_TOOLS = {
    "render_structure_svg": render_structure_svg,
    "describe_selection": describe_selection_tool,
    "render_reactant_scheme": render_reactant_scheme,
}

PAGE = """<!doctype html>
<html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>chemcard mock host</title>
<style>
  :root {{ --bg:#fafaf9; --fg:#1c1c1c; --log:#666; color-scheme: light; }}
  :root[data-theme="dark"] {{ --bg:#1f1f1e; --fg:#ececec; --log:#999; color-scheme: dark; }}
  body {{ margin:0 auto; max-width:720px; padding:24px 16px; background:var(--bg); color:var(--fg);
         font:16px/1.6 system-ui, sans-serif; }}
  iframe {{ display:block; width:100%; border:0; height:0; background:transparent; }}
  .log {{ font:11px/1.4 ui-monospace, monospace; color:var(--log); white-space:pre-wrap; margin:0 0 12px; }}
  #composer {{ position:sticky; bottom:0; background:var(--bg); border-top:1px solid #8884; padding:8px 0;
              font:12px ui-monospace, monospace; white-space:pre-wrap; }}
  button {{ font:inherit; }}
  body.demo .log, body.demo p:has(#theme), body.demo #composer {{ display:none; }}
  body.demo p {{ margin: 14px 0 6px; }}
</style></head><body class="{body_class}">
<p><button id="theme">切换亮/暗主题</button></p>
{blocks}
<div id="composer">输入框（ui/message 会出现在这里）：</div>
<script>
const WIDGET = {widget};
const RESULTS = {results};
const INPUTS = {inputs};
const MODES = {modes};
const CAPS = {{ serverTools: {{}}, message: {{ text: {{}} }}, updateModelContext: {{ text: {{}} }} }};
let theme = "light";
const frames = [];
document.querySelectorAll("iframe").forEach((frame, i) => {{
  const log = document.getElementById("log" + i);
  const say = (s) => log.textContent += s + "\\n";
  const post = (m) => frame.contentWindow.postMessage(m, "*");
  frames.push(frame);
  window.addEventListener("message", async (ev) => {{
    if (ev.source !== frame.contentWindow) return;
    const m = ev.data;
    if (m.method === "ui/notifications/size-changed") {{ frame.style.height = m.params.height + "px"; return; }}
    say("← " + m.method);
    if (m.method === "ui/initialize") {{
      post({{ jsonrpc:"2.0", id:m.id, result:{{ protocolVersion:"2026-01-26", hostInfo:{{name:"mock-host",version:"0"}},
        hostCapabilities: CAPS, hostContext:{{ theme, displayMode:"inline" }} }} }});
    }} else if (m.method === "ui/notifications/initialized") {{
      if (MODES[i] === "silent") return;
      post({{ jsonrpc:"2.0", method:"ui/notifications/tool-input", params:{{ arguments: INPUTS[i] }} }});
      const params = Object.assign({{}}, RESULTS[i]);
      if (MODES[i] === "strip_meta") delete params._meta;
      post({{ jsonrpc:"2.0", method:"ui/notifications/tool-result", params }});
      say("→ tool-input, tool-result" + (MODES[i] === "strip_meta" ? "（已去掉 _meta）" : ""));
    }} else if (m.method === "tools/call") {{
      const r = await fetch("/call", {{ method:"POST", body: JSON.stringify(m.params) }}).then(r => r.json());
      post({{ jsonrpc:"2.0", id:m.id, result: r }});
      say("→ tools/call " + m.params.name + (m.params.name === "describe_selection" ? " " + JSON.stringify(r.structuredContent) : ""));
    }} else if (m.method === "ui/update-model-context") {{
      post({{ jsonrpc:"2.0", id:m.id, result:{{}} }});
      say("   上下文 = " + m.params.content.map(c => c.text).join(" "));
    }} else if (m.method === "ui/message") {{
      post({{ jsonrpc:"2.0", id:m.id, result:{{}} }});
      document.getElementById("composer").textContent = "输入框（ui/message 会出现在这里）：\\n" + m.params.content.map(c => c.text).join("");
    }} else if (m.id !== undefined) {{
      post({{ jsonrpc:"2.0", id:m.id, result:{{}} }});
    }}
  }});
  frame.srcdoc = WIDGET;
}});
document.getElementById("theme").onclick = () => {{
  theme = theme === "light" ? "dark" : "light";
  document.documentElement.dataset.theme = theme;
  frames.forEach(f => f.contentWindow.postMessage({{ jsonrpc:"2.0", method:"ui/notifications/host-context-changed", params:{{ theme }} }}, "*"));
}};
</script></body></html>
"""


def _dump(result) -> dict:
    return result.model_dump(by_alias=True, exclude_none=True, mode="json")


def _js(value: object) -> str:
    """JSON safe to inline in <script> (the widget itself contains '</script>')."""
    return json.dumps(value).replace("</", "<\\/")


async def _read_widget() -> str:
    contents = list(await mcp.read_resource(CARD_URI))
    return contents[0].content


def build_page(demo: bool = False) -> str:
    widget = anyio.run(_read_widget)
    blocks, results, inputs, modes = [], [], [], []
    for i, (text, tool, args, mode) in enumerate(DEMO if demo else SAMPLES):
        results.append(_dump(DRAW_TOOLS[tool](**args)))
        inputs.append(json.loads(json.dumps(args, default=lambda o: o.model_dump(exclude_none=True))))
        modes.append(mode)
        blocks.append(
            f"<p>{html.escape(text)}</p>\n"
            f'<iframe sandbox="allow-scripts" title="card {i}"></iframe>\n'
            f'<div class="log" id="log{i}"></div>'
        )
    return PAGE.format(blocks="\n".join(blocks), widget=_js(widget), results=_js(results),
                       inputs=_js(inputs), modes=_js(modes), body_class="demo" if demo else "")


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        demo = "demo=1" in self.path
        body = build_page(demo).encode()  # rebuilt per load, so edits show up on refresh
        self._send(200, "text/html; charset=utf-8", body)

    def do_POST(self):
        req = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
        tool = APP_TOOLS.get(req.get("name"))
        result = _dump(tool(**req.get("arguments", {}))) if tool else {"isError": True}
        self._send(200, "application/json", json.dumps(result).encode())

    def _send(self, code, ctype, body):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"mock host: http://127.0.0.1:{args.port}/", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
