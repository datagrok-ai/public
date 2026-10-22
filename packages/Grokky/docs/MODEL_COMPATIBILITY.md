# Grokky model compatibility

Which model families and providers Grokky works on, and the evidence for each claim. Measured 2026-10-22 on
branch `claude/grokky-providers`, with real model calls from the AI panel in a browser on a local stand.

## Summary

Grokky runs on Claude and GPT. Through Anthropic and OpenAI, the cheapest model of each family passes the whole
AI-panel BDD suite, 22 scenarios, on the final runtime image: computed answers, a dashboard with a follow-up, a
viewer refined over three turns, filtering and extracting rows, a pivot table checked cell by cell, summary tables
in their own views, server work (a project saved, a folder listed, a delete that waits for approval), docs
answers, asking the user, sandbox Python, a stopped job, a reload, a lost agent process, a `!command` and a
security probe. OpenRouter passes five core scenarios on both families. Bedrock passed its live checks but
its browser runs are blocked on a key; Azure and Databricks are checked against stubs and, for Databricks, the
route and sign-in only. The runtime keeps no model list: other families (GPT OSS, Llama, Qwen) can be named
through OpenRouter or a custom endpoint, but were not tested.

Making GPT work well took tuning, not code per provider: it read the js-api source where a skill had the answer,
explored the view before answering a plain question, and computed statistics over missing values. The prompt and
skill changes below fixed that for both families; over the whole suite GPT's source reads fell from 24 to 6 and
Claude's from 8 to 0.

## Status

| Family | Route | Status | Evidence |
|---|---|---|---|
| Claude | Anthropic | **Validated** | Full BDD suite 22/22 on Haiku 5.5 (final image, and the image before it); 7/7 core in four earlier runs; 18/18 on Sonnet 5.5 before the suite grew; live checks 5/5 |
| Claude | OpenRouter | **Validated** (5 core scenarios) | BDD 5/5 after the shared fixes; live checks 5/5 |
| Claude | Custom endpoint, `anthropic-messages` | Partially validated | BDD 2/2 (computed answer, dashboard); live checks 5/5 |
| Claude | Amazon Bedrock | Partially validated | Live checks 5/5; BDD blocked: the short-term key ended with its console session |
| Claude | Microsoft Azure | Not tested | Stubs under Azure's host names only |
| Claude | Databricks | Not tested | Route and sign-in only, on a Free Edition workspace without Claude |
| GPT | OpenAI | **Validated** | Full BDD suite 22/22 on GPT-6 Luna (final image, and the image before it); live checks 5/5 |
| GPT | OpenRouter | **Validated** (5 core scenarios) | BDD 5/5 after the shared fixes; live checks 5/5 |
| GPT | Custom endpoint, `openai-responses` | Partially validated | BDD 2/2 after the shared fixes; live checks 5/5 |
| GPT | Custom endpoint, `openai-chat-completions` | Partially validated | Live only: tools, bash and streaming pass; GPT-6 tool calls with reasoning are refused by OpenAI, and structured output fails on gateways that refuse an output format beside tools |
| GPT | Amazon Bedrock | Partially validated | Live checks 5/5, structured output included; BDD blocked as for Claude |
| GPT | Microsoft Azure | Not tested | Stubs only |
| GPT | Databricks | Not tested | Route and sign-in only; Free Edition hosts no GPT-6 |
| Open-weight (GPT OSS, Llama, Qwen) | OpenRouter, custom endpoint | Not tested | Can be named in `model`; called in the Responses or the endpoint's format |

**Validated** means the BDD scenarios that apply passed in a recorded browser run, tools, skills and multi-turn
tasks included, on top of the live checks. Each family was tested on its cheapest model only (Claude Haiku 5.5,
GPT-6 Luna); larger models of a family are not covered by these runs.

## Configuration

| Route | Settings (`provider`, `model`) |
|---|---|
| Anthropic | default provider, `claude-haiku-5-5` |
| OpenAI | `OpenAI`, `gpt-6-luna` |
| OpenRouter | `OpenRouter`, `anthropic/claude-haiku-5.5` or `openai/gpt-6-luna` |
| Bedrock | `Bedrock`, region `us-east-2`, `global.anthropic.claude-haiku-5-5` or `global.openai.gpt-6-luna` |
| Custom endpoint | `Custom`, `anthropic-messages` at `https://api.anthropic.com` or `openai-responses` at `https://api.openai.com/v1` |

The stand ran this branch's plugin and runtime image (rebuilt after each change), Datagrok core at the commit its
worktree pinned, one browser worker, and video at 1080p. Before each run the runtime itself was asked which
models it served (its `models` reply), and every video carries that model in its caption. The runs before the
full suite waited about a minute per turn for the stand's admission queue, which had no worker; the chat function
no longer declares a queue, so turns now go straight to the runtime.

## Scenarios

The whole suite is the 18 scenarios of `bdd/features/ai-panel/`: the seven core scenarios below and the rest
(server work, session, security, two known platform issues). The 22-scenario runs also had four exploratory
analysis scenarios, run locally and not kept in the repo:

| Scenario | What it shows |
|---|---|
| A question about the table is answered with the value computed from it | a tool computes the answer |
| A small dashboard is built from one request, then changed by a follow-up | multi-step tool use, context across turns, the result in the UI |
| A how-to question is answered from the documentation, with a link to it | the docs index, a help page and its link |
| The table view's own commands color-code a column and hide another | finding and calling the view's functions |
| A request for data the table does not have changes nothing | handling a request that cannot be met |
| The AI asks which chart to add and adds the one picked | the ask-the-user round trip |
| Python in the AI's sandbox computes over a file from the server | skills, the MCP file tools and the sandbox |
| *Exploratory:* A histogram is refined over three turns, every change landing on the same viewer | one viewer changed turn after turn, never duplicated |
| *Exploratory:* Rows are filtered, a subgroup selected, and the selection extracted into a table of its own | filter, select and extract, checked at 225 rows |
| *Exploratory:* A pivot table gives the mean weight of every race, one column per sex | every cell against an independent aggregation |
| *Exploratory:* A per-race summary is charted in its own view | a new table, its view and a chart in it |

The table view's own prompt handler takes a request that names a single viewer, and a viewer's quick commands
take short edits ("split it by SEX", "20 bins"), before any model sees them; a prompt meant for the model must be
worded so the platform passes it on. To repeat a run:

```bash
GROKKY_BDD_LIVE=1 DATAGROK_SERVER=<alias> DATAGROK_URL=<stand> npx grok-bdd run --workers 1 -g "<scenario title>"
```

The videos and the per-turn tool logs of these runs came from local changes to the BDD harness, not in the repo.

## Tuning: skills and tools instead of the source

The whole suite, from the runtime's own events (one run each):

| Run | Passed | Source reads and searches | Skill loads | Platform tools | Model calls | Model API time |
|---|---|---|---|---|---|---|
| GPT-6 Luna, first full run | 15/22 | 18 | 12 | 90 | 117 | 721 s |
| GPT-6 Luna, after the fixes | 22/22 | 24 | 16 | 133 | 205 | 718 s |
| GPT-6 Luna, final: + skill routing | **22/22** | **6** | 19 | 115 | 165 | 496 s |
| Claude Haiku 5.5, first full run | 17/22 | 5 | 10 | 60 | 99 | 316 s |
| Claude Haiku 5.5, after the fixes | 22/22 | 8 | 11 | 86 | 133 | 401 s |
| Claude Haiku 5.5, final: + skill routing | **22/22** | **0** | 14 | 76 | 121 | 336 s |

The first run's failures that were the model's (the rest were prompts the platform answered itself and a known
issue tagged on the wrong line) and what fixed them, the same for every model:

- **Statistics over missing values** (`datagrok-exec`, inlined in every prompt): GPT computed a median from the
  raw array, where a missing number is a placeholder, and answered 165.8 instead of 168.5. A number from a
  column now comes from `col.stats`, which skips missing values.
- **Per-group summaries** (`datagrok-df-and-columns`): no skill named `groupBy`, so GPT searched the source until
  the turn timed out. The skill now gives the recipe, and its description names statistics and per-group summary
  tables, which is what makes a model open it.
- **Box plot and pivot table properties** (`datagrok-viewers`): the skill gave the box plot a
  `categoryColumnName` it does not have, so a model that followed it got the default split (DIS_POP for RACE) and
  reported success; it named no pivot table property at all. Both now list the real ones, and the skill says how to
  read a viewer's current setting back (`viewer.props.<name>`).

Earlier, on the seven core scenarios:

| Run | Passed | Source reads and searches | Skill loads | Platform tools | Model calls | Model API time |
|---|---|---|---|---|---|---|
| GPT-6 Luna, before | 7/7 | 6 | 9 | 38 | 60 | 240 s |
| GPT-6 Luna, lookup order | 7/7 | 1 | 8 | 37 | 61 | 256 s |
| GPT-6 Luna, + read-only questions | 6/7 | 1 | 6 | 32 | 44 | 210 s |
| GPT-6 Luna, + missing values | 7/7 | **0** | 7 | 43 | 60 | 277 s |
| Claude Haiku 5.5, before | 7/7 | 0 | 5 | 21 | 36 | 125 s |
| Claude Haiku 5.5, after all three | 7/7 | 0 | 6 | 26 | 39 | 134 s |

- **Lookup order** (`prompts.ts`, "Don't invent names"): names come from skills, then `skill_view`, the view's
  functions and MCP discovery; a skill's name and arguments are final; the source is a last resort, one file, no
  search.
- **Read-only questions** (`prompts.ts`, "View functions"): the view's functions are for requests that change the
  view; a question about the data is one read-only `datagrok_exec`. GPT's table question went from about ten tool
  calls to two.
- **Missing values in a filter** (`datagrok-filtering`): `t.rows.removeWhereIdx((i) => t.getCol('X').isNone(i))`.

GPT still runs a multi-step task as more, smaller `datagrok_exec` calls than Claude, which costs time, not
correctness; that needs the benchmark (`docs/BENCHMARK.md`) with repeated runs before it is tuned.

## Defects found and fixed

| Found | Fix |
|---|---|
| A delete's confirm token never expired, so Claude Haiku reused one the user had declined and deleted a project with no card | the runtime spends a token once, and only after the user approved its card in the same turn |
| The skills' box plot property did not exist; no pivot table properties; no way to read a setting back | `datagrok-viewers` |
| No skill covered per-group aggregation or said statistics must skip missing values | `datagrok-df-and-columns`, `datagrok-exec` |
| Every turn waited a minute for an admission queue with no worker | the chat function no longer declares a queue |
| A plain-text link in an answer was not clickable (Claude Haiku wrote one) | the panel turns bare URLs into links (`linkBareUrls`) |
| Nothing told the model how a help page maps to its URL; GPT built `datagrok.ai/help.datagrok.ai/…` | the help index states the mapping |
| GPT read the js-api source, explored the view for plain questions, missed missing values | the tuning changes above |
| The stand served a js-api bundle older than its Dart client, so `addNewCalculated` failed for any model | stand only: the bundle rebuilt (`pnpm turbo run build --filter=datagrok-api --force`) |

Found in the platform and not fixed here (core): the table view's prompt handler takes a request that names one
viewer even when it asks for more, adds that viewer and reports the request handled ("make a summary table, open
it and chart it as bars" adds a bar chart to the current view, with no table); and "put a histogram of AGE on this
view, and tell me…" adds a bar chart.

Documented in `CLAUDE_CODE_FLOW.md`: OpenRouter needs purchased credits, not only a key limit; a Bedrock
short-term key ends with its console session; OpenAI's own Chat Completions refuses GPT-6 tool calls with a
reasoning effort, so GPT-6 runs in `openai-responses`; a gateway may refuse an output format beside tools.

## Cost and time

Model cost from the runtime's metrics, over 304 scenario runs: Anthropic $0.57, OpenAI $0.67, OpenRouter $0.09,
custom endpoints $0.03, Bedrock $0 (every call refused). A full suite of 22 scenarios on the final image cost
$0.094 and 15 minutes on GPT-6 Luna, $0.076 and 12 minutes on Claude Haiku 5.5.

## Gaps and next tests

1. Bedrock in the browser with a long-term API key.
2. Azure and Databricks with live models (a Foundry resource; a workspace that serves the `system.ai` Claude and
   GPT models).
3. The full suite on OpenRouter and the custom endpoint, not only the core scenarios.
4. The larger models of each family (Claude Sonnet and Opus 5.5, GPT-6 Sol and 6.1 Sol) in the browser.
5. Repeated runs per scenario, to tell variance from a defect.
6. GPT's last source reads: how to open a new table view and read its table.
7. The 52-prompt benchmark per model, to tune GPT's step count with numbers behind it.
