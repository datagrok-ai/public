# Grokky model compatibility

Which model families and providers Grokky works on, and the evidence for each claim. Measured 2026-10-22 on
branch `claude/grokky-providers`, with real model calls from the AI panel in a browser on a local stand.

## Summary

Grokky runs on Claude and GPT. Through Anthropic, OpenAI and OpenRouter, the cheapest model of each family
passes every representative BDD scenario: a computed answer, a multi-step dashboard with a follow-up, a docs
answer with its link, the view's own commands, a request for missing data, asking the user, and sandbox Python
over a server file. Bedrock passed its live checks but its browser runs are blocked on a key; Azure and
Databricks are checked against stubs and, for Databricks, the route and sign-in only. Only Claude and GPT are in
Grokky's model catalog; open-weight families (GPT OSS, Llama, Qwen) are not supported.

GPT needed tuning, not code changes per provider: it read the js-api source instead of using Grokky's skills,
explored the view before answering a plain question, and missed how Datagrok marks a missing value. Three prompt
and skill changes fixed that (see [Tuning](#tuning-skills-and-tools-instead-of-the-source)); both families pass
everything they ran, with no source reads.

## Status

| Family | Route | Status | Evidence |
|---|---|---|---|
| Claude | Anthropic | **Validated** | BDD 7/7 on Haiku 5.5 (four runs); full suite 18/18 on Sonnet 5.5; live checks 5/5 |
| Claude | OpenRouter | **Validated** (5 core scenarios) | BDD 5/5 after the shared fixes; live checks 5/5 |
| Claude | Custom endpoint, `anthropic-messages` | Partially validated | BDD 2/2 (computed answer, dashboard); live checks 5/5 |
| Claude | Amazon Bedrock | Partially validated | Live checks 5/5; BDD blocked: the short-term key ended with its console session |
| Claude | Microsoft Azure | Not tested | Stubs under Azure's host names only |
| Claude | Databricks | Not tested | Route and sign-in only, on a Free Edition workspace without Claude |
| GPT | OpenAI | **Validated** | BDD 7/7 on GPT-6 Luna in four of five runs (the fifth missed missing values, since fixed); live checks 5/5 and 19/20 of the wider suite |
| GPT | OpenRouter | **Validated** (5 core scenarios) | BDD 5/5 after the shared fixes; live checks 5/5 |
| GPT | Custom endpoint, `openai-responses` | Partially validated | BDD 2/2 after the shared fixes; live checks 5/5 |
| GPT | Custom endpoint, `openai-chat-completions` | Partially validated | Live only: tools, bash and streaming pass; GPT-6 tool calls with reasoning are refused by OpenAI, and structured output fails on gateways that refuse an output format beside tools |
| GPT | Amazon Bedrock | Partially validated | Live checks 5/5, structured output included; BDD blocked as for Claude |
| GPT | Microsoft Azure | Not tested | Stubs only |
| GPT | Databricks | Not tested | Route and sign-in only; Free Edition hosts no GPT-6 |
| Open-weight (GPT OSS, Llama, Qwen) | any | Not supported | Not in the model catalog |

**Validated** means the BDD scenarios that apply passed in a recorded browser run, tools, skills and a
multi-turn task included, on top of the live checks. Each family was tested on its cheapest model only (Claude
Haiku 5.5, GPT-6 Luna); larger models of a family are not covered by these runs.

## Configuration

| Route | Settings (`provider`, `model`, `fastModel`) |
|---|---|
| Anthropic | default provider, `claude-haiku-5-5` / `claude-haiku-5-5` |
| OpenAI | `OpenAI`, `gpt-6-luna` / `gpt-6-luna` |
| OpenRouter | `OpenRouter`, `anthropic/claude-haiku-5.5` or `openai/gpt-6-luna`, the same as fast |
| Bedrock | `Bedrock`, region `us-east-2`, `global.anthropic.claude-haiku-5-5` or `global.openai.gpt-6-luna` |
| Custom endpoint | `Custom`, `anthropic-messages` at `https://api.anthropic.com` or `openai-responses` at `https://api.openai.com/v1` |

The stand ran this branch's plugin and runtime image (rebuilt after each change), Datagrok core at the commit its
worktree pinned, one browser worker, and video at 1080p. Before each run the runtime itself was asked which
models it served (its `models` reply), and every video carries that model in its caption. Every turn waited about
a minute for the stand's admission queue, which had no worker, so wall-clock times overstate the model's own
time; the model time below comes from each turn's metrics.

## Scenarios

The seven scenarios are the existing ones in `bdd/features/ai-panel/`, reused unchanged:

| Scenario | What it shows |
|---|---|
| A question about the table is answered with the value computed from it | a tool computes the answer |
| A small dashboard is built from one request, then changed by a follow-up | multi-step tool use, context across turns, the result in the UI |
| A how-to question is answered from the documentation, with a link to it | the docs index, a help page and its link |
| The table view's own commands color-code a column and hide another | finding and calling the view's functions |
| A request for data the table does not have changes nothing | handling a request that cannot be met |
| The AI asks which chart to add and adds the one picked | the ask-the-user round trip |
| Python in the AI's sandbox computes over a file from the server | skills, the MCP file tools and the sandbox |

No scenario was added: structured output, streaming, recovery after a lost agent process and the bad-key errors
are covered by the runtime's live suite on each route, and the full BDD suite (18 scenarios) ran on Claude
Sonnet 5.5. To repeat a run:

```bash
GROKKY_BDD_LIVE=1 DATAGROK_SERVER=<alias> DATAGROK_URL=<stand> GROKKY_BDD_TURNS=<file.jsonl> \
  npx grok-bdd run --video on --workers 1 -g "<scenario title>"
```

`--video on` films each run into `bdd/test-results/videos/`; `GROKKY_BDD_TURNS` logs each turn's model, tool
calls, answer, cost and tokens.

## Tuning: skills and tools instead of the source

Tool calls over the seven scenarios, from the runtime's own events (one run each; a single run varies by about
±10 calls on the dashboard):

| Run | Passed | Source reads and searches | Skill loads | Platform tools | Model calls | Model API time |
|---|---|---|---|---|---|---|
| GPT-6 Luna, before | 7/7 | 6 | 9 | 38 | 60 | 240 s |
| GPT-6 Luna, lookup order | 7/7 | 1 | 8 | 37 | 61 | 256 s |
| GPT-6 Luna, + read-only questions | 6/7 | 1 | 6 | 32 | 44 | 210 s |
| GPT-6 Luna, + missing values | 7/7 | **0** | 7 | 43 | 60 | 277 s |
| Claude Haiku 5.5, before | 7/7 | 0 | 5 | 21 | 36 | 125 s |
| Claude Haiku 5.5, after all three | 7/7 | 0 | 6 | 26 | 39 | 134 s |

- **Lookup order** (`prompts.ts`, "Don't invent names"): the rule sent models to `workspace/js-api/src/`
  whenever a name was not in an inlined skill. GPT searched it for `addNewCalculated` and `class Viewer` on the
  dashboard and called the method with the source's full overload. Names now come from skills, then
  `skill_view`, the view's functions and MCP discovery; a skill's name and arguments are final; the source is a
  last resort, one file, no search.
- **Read-only questions** (`prompts.ts`, "View functions"): "FIRST call `list_view_functions`" made GPT explore
  the view, call its functions and verify before answering which race weighs most; now that rule applies to
  requests that change the view, and a question about the data is one read-only `datagrok_exec`. GPT's table
  question went from about ten tool calls to two.
- **Missing values** (`datagrok-filtering` skill): GPT tested missing heights with its own check and reported
  "none were missing" (128 were). The skill now gives `t.rows.removeWhereIdx((i) => t.getCol('X').isNone(i))`.

GPT still runs a multi-step task as more, smaller `datagrok_exec` calls than Claude (12 against 7 on the
dashboard's first turn). That costs time, not correctness, and needs the benchmark (`docs/BENCHMARK.md`) with
repeated runs before it is tuned.

## Defects found and fixed

| Found | Fix |
|---|---|
| `grok-bdd run --video on` was refused by Playwright and could not film the shared page | `libraries/bdd`: the flag passes `BDD_VIDEO`, the worker's page is filmed and written at feature end |
| A plain-text link in an answer was not clickable (Claude Haiku wrote one) | the panel turns bare URLs into links (`linkBareUrls`) |
| Nothing told the model how a help page maps to its URL; GPT built `datagrok.ai/help.datagrok.ai/…` | the help index states the mapping |
| GPT read the js-api source, explored the view for plain questions, missed missing values | the three tuning changes above |
| The stand served a js-api bundle older than its Dart client, so `addNewCalculated` failed for any model | stand only: the bundle rebuilt (`pnpm turbo run build --filter=datagrok-api --force`) |
| Two live-suite checks assumed Claude's behaviour (isolation, schema patterns) | the local suite judges both families |

Documented in `CLAUDE_CODE_FLOW.md`: OpenRouter needs purchased credits, not only a key limit; a Bedrock
short-term key ends with its console session; OpenAI's own Chat Completions refuses GPT-6 tool calls with a
reasoning effort, so GPT-6 runs in `openai-responses`; a gateway may refuse an output format beside tools.

## Cost and time

Model cost from the runtime's metrics, over 105 scenario runs: Anthropic $0.20, OpenAI $0.23, OpenRouter $0.09,
custom endpoints $0.04, Bedrock $0 (every call refused). On the seven scenarios before tuning, GPT-6 Luna cost
$0.043 against Claude Haiku's $0.031, with about twice the model calls and API time.

## Gaps and next tests

1. Bedrock in the browser with a long-term API key.
2. Azure and Databricks with live models (a Foundry resource; a workspace that serves the `system.ai` Claude and
   GPT models).
3. The larger models of each family (Claude Sonnet and Opus 5.5, GPT-6 Sol and 6.1 Sol) in the browser.
4. Repeated runs per scenario, to tell variance from a defect.
5. The 52-prompt benchmark per model, to tune GPT's step count with numbers behind it.
6. The custom endpoint's docs and view-command scenarios, and a Chat Completions gateway in the browser.
