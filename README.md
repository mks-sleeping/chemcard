# chemcard

> Inline 2D chemical structures for AI chat explanations — an MCP App (RDKit + MCP Apps) that lets Claude or ChatGPT drop clean structure figures, reaction steps and highlights between paragraphs, like equations in a textbook.

让 AI 在讲解化学时，在文字之间插入干净、正确、紧凑的 2D 结构图：像讲义里的公式一样，而不是嵌入一个完整的化学编辑器。

![chemcard 效果：讲解与结构图交替，高亮正在讲的基团](docs/images/demo.png)

> **早期实验项目（v0.3）。** 已在 Claude 桌面版中实测；ChatGPT 接入方式已按文档实现，但尚未实测。

## 能做什么

- **画分子**（`draw_molecule`）：模型给 SMILES，服务器用 RDKit 画图。键长固定，小分子画得小，大分子不会被压缩；桥环（如托品烷）使用环系模板，画法清晰。
- **画反应步骤**（`draw_reaction_step`）：一步一张「A + B → C」，条件写在箭头上方，可逆反应用 ⇌。讲机理时一步一张，图下紧跟这一步的讲解。
- **高亮**：用 SMARTS（`"C=O"`）或原子映射号（`[NH2:1]`）给原子和键加上彩色底色，标出正在讲的部分。
- **点选原子**：在卡片上点选原子，自动转成片段 SMILES 和带编号的 SMILES，推给模型或填入输入框；模型可以用同样的编号回过头来高亮。
- **检查机理**（`check_mechanism`，不显示卡片）：检查每个物种的价态、意外的自由基（例如 `C[N+](C)C` 漏写了 H）、不稳定的碳正离子、酸碱条件矛盾，以及每一步的原子与电荷是否守恒、前后步骤是否衔接。明显的笔误会被直接拒绝，模型收到原因后自行修正。

## 快速开始

需要 [uv](https://docs.astral.sh/uv/) 和 Python ≥ 3.11。

```bash
git clone <this repo> chemcard && cd chemcard
uv sync
uv run pytest
```

### 接入 Claude 桌面版（本地运行，无需部署）

在 `~/Library/Application Support/Claude/claude_desktop_config.json` 的 `mcpServers` 中加入下面的配置。`uv` 请写**绝对路径**（可用 `which uv` 查看），因为从 Dock 启动的应用读不到终端里的 PATH。

```json
"chemcard": {
  "command": "/absolute/path/to/uv",
  "args": ["--directory", "/absolute/path/to/chemcard", "run", "chemcard"]
}
```

完全退出 Claude（⌘Q）再打开，新开对话问「画出阿司匹林，标出会水解的部分」。

**想让模型主动画图？** Claude 的连接器工具是按需加载的；问题本身没有提到「画」或「结构」时，模型往往不会去加载 chemcard。建议在 设置 → 个人资料 → 个人偏好（或项目说明）里加一句：

```text
讲解化学时，只要涉及具体的分子、官能团或反应机理，就用 chemcard 连接器把结构画出来（draw_molecule / draw_reaction_step），并用高亮标出正在讲的部分。
```

### 接入 ChatGPT（需要公网 HTTPS，未实测）

```bash
CHEMCARD_ALLOWED_HOSTS=xxxx.ngrok-free.app uv run chemcard --http   # 监听 http://127.0.0.1:8000/mcp
ngrok http 8000
```

`CHEMCARD_ALLOWED_HOSTS` 用于放行隧道域名，否则 MCP SDK 的 DNS rebinding 防护会返回 `421 Invalid Host header`。然后在 ChatGPT 中选择「Add custom MCP server」，填入 `https://xxxx.ngrok-free.app/mcp`。

### 本地预览（不需要任何宿主）

```bash
uv run python scripts/mock_host.py
```

打开 http://127.0.0.1:8765/ ：这是一个模拟聊天页，用真实卡片跑一遍 MCP Apps 握手，可以切换亮/暗主题、点选原子，每张卡片收发的消息都会显示出来。打开 http://127.0.0.1:8765/?demo=1 是 README 截图用的干净版本。

## 工具一览

| 工具 | 谁调用 | 作用 |
|---|---|---|
| `draw_molecule` | 模型 | 画一个分子，可带 `highlights` |
| `draw_reaction_step` | 模型 | 画一步「A + B → C」，`highlights` 用 `side` 和 `index` 指明标在哪个物种上 |
| `check_mechanism` | 模型 | 画多步机理之前整体检查一遍，不显示卡片 |
| `render_structure_svg` | 仅卡片 | 宿主没有把图交给卡片时，卡片自己取图 |
| `describe_selection` | 仅卡片 | 把点选的原子转成片段 SMILES 和带编号的 SMILES |

## 已知限制

- **模型不一定会主动画图**：需要用户侧的提示，比如个人偏好、项目说明或 Skill，见上文。
- **化学检查只查账面错误**：查得出价态、守恒和前后衔接的问题，判断不了「这种条件下会不会这样反应」，也判断不了选择性。
- **还没有弯箭头**（电子推动箭头），在计划中。
- **Claude 是否支持 `ui/update-model-context` 未确认**：如果不支持，点选原子后需要点「提问」，模型才知道你选了什么。
- **ChatGPT 未实测。**
- **重新打开历史对话时，chemcard 服务器需要在运行**：Claude 不把 `_meta` 交给卡片，卡片要回服务器取图。

## 踩坑记录

如果你也在做 MCP Apps，这些可能帮你省点时间（均为 2026-11 在 Claude 桌面版中观察到的现象）：

- **Claude 会丢掉工具参数结构里的 `$ref`**：嵌套的模型字段会变成 `items: {}`，模型看不到这些字段，自然也不会填。解决办法是注册工具后把参数结构展开（见 `server.flatten_schema`），并让每个工具的嵌套深度不超过一个已确认能完整显示的工具。
- **Claude 不把工具结果里的 `_meta` 交给卡片**：卡片拿不到 SVG，就会显示为空白。解决办法是提供一个只给卡片调用的工具（`visibility: ["app"]`），让卡片通过 `tools/call` 自己取图。
- **卡片向宿主报告高度 0 会导致空白卡片**：内容出现之前不要发送 `size-changed`。
- **RDKit 的自适应画布（`MolDraw2DSVG(-1, -1)`）会忽略 `fixedBondLength`**：先画一遍测出实际键长，再按比例缩放坐标后重画。
- **CoordGen 给单独一个 `[H+]` 算出的坐标是 NaN**，画出来是 0×0 的空图；桥环骨架用 CoordGen 也画得很别扭。这两种情况改用 RDKit 自带的布局算法并开启 `useRingTemplates`。
- **RDKit 的 SVG 声明编码为 `iso-8859-1`**：图例里的中文会变成乱码，所以图注改放在 HTML 里。
- **RDKit 会把 `C[N+](C)C` 当成自由基阳离子，而不会报错**：如果要防住模型漏写 H，需要自己加规则检查。
- **模型每次调用后立刻读到的返回文字，比工具描述更能左右它的行为**：「先讲这一步，再画下一步」这句提醒，写在返回文字里比写在描述里有效得多。

## 架构

```
src/chemcard/
  chem/            ① 纯化学核心：不 import mcp，可单独测试
    render.py         SMILES → SVG，分子式，标准化 SMILES，原子坐标
    highlight.py      高亮解析（SMARTS / 原子映射号 → 标准原子序号）与选区描述
    reaction.py       反应步骤的原子与电荷守恒检查
    lint.py           化学规则检查
    mechanism.py      整条机理检查（规则检查 + 守恒 + 前后衔接）
  server.py        ② MCP 适配层：注册工具和资源、组装返回值，不做化学计算
  descriptions.py     模型读到的所有文字（工具描述、服务器说明、调用后的提醒）
  widgets/card.html ③ 卡片：原生 HTML/CSS/JS，只负责显示，不解析化学
tests/             pytest（100+ 个用例）
evals/             触发评测：哪些问题应该画图、图文是否穿插，以及实测记录
scripts/mock_host.py  本地模拟宿主
```

工具返回值分三部分，各自给不同的读者：

| 字段 | 谁能看到 | 放什么 |
|---|---|---|
| `content` | 模型和卡片（不支持 UI 的客户端只显示它） | 一句文字说明，以及给模型的提醒 |
| `structuredContent` | 模型和卡片 | 精简信息：标准化 SMILES、分子式、图注、高亮 |
| `_meta` | 只有卡片 | SVG 本体和原子坐标，不占用模型上下文 |

## 路线

| 版本 | 内容 |
|---|---|
| v0.1 ✓ | 单分子卡片 |
| v0.2 ✓ | 反应步骤卡片，以及机理检查（`check_mechanism`） |
| v0.3 ✓ | 高亮；点选原子回传给模型 |
| v0.4 | 弯箭头：在反应物上叠加电子推动箭头，用原子映射号指定起点和终点 |
| v0.5 | 电子推动校验：按箭头推算产物，并与给出的产物比对 |

## 致谢与许可

本项目使用 [MIT 许可证](LICENSE)。

- [RDKit](https://www.rdkit.org/)（BSD-3-Clause）：所有化学解析和绘图。
- [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk) 与 [MCP Apps 规范](https://github.com/modelcontextprotocol/ext-apps)。
- 思路参考（没有复制任何代码）：[ChemCP](https://github.com/scottmreed/ChemCP)、[chemdraw-mcp](https://github.com/jurimaxam-dotcom/chemdraw-mcp)、[ChemKit](https://github.com/Oli-CodeHub/ChemKit)，以及 OpenAI 的 [Apps SDK 示例](https://github.com/openai/openai-apps-sdk-examples)。

本项目与 RDKit、Anthropic、OpenAI、Revvity（ChemDraw）均无关联。
