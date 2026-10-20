# Claude Code Flow

This document describes how the Claude agent runtime is deployed and integrated into Datagrok: how the server-side
containers are set up, how messages travel between the browser and the containers, and how Claude's
responses are turned into live Datagrok UI.

---

## Architecture Overview

Three components work together:

```
Browser (TypeScript)
  └── ClaudeRuntimeClient  ──WebSocket──▶  claude-runtime container  ──MCP──▶  mcp-server container
                                                  │                                    │
                                     one agent process per session               grok.dapi.*
                                  (Hermes Agent, JSON lines on stdio)     (Datagrok API over HTTP)
```

| Component                  | Location                                                          | Port |
|----------------------------|-------------------------------------------------------------------|------|
| Browser client             | [`src/claude/runtime-client.ts`](../src/claude/runtime-client.ts) | —    |
| `claude-runtime` container | [`dockerfiles/claude-runtime/`](../dockerfiles/claude-runtime/)   | 5355 |
| `mcp-server` container     | [`dockerfiles/mcp-server/`](../dockerfiles/mcp-server/)           | 3003 |

---

## Deployment

### Containers

Both containers live under [`dockerfiles/`](../dockerfiles/) and are managed by the Datagrok platform's
Docker subsystem. The platform discovers and starts them automatically when the Grokky package is published.

**`claude-runtime`** ([`Dockerfile`](../dockerfiles/claude-runtime/Dockerfile)):

- Base image: `node:22-trixie-slim`, pinned by digest.
- Clones the public Datagrok repo into `/workspace` (and the help repo into `/workspace/help`) at build
  time — this gives Claude read access to the JS API source, packages, samples, and documentation
  without fetching them at runtime.
- Runs `dist/server.js` (compiled from [`src/server.ts`](../dockerfiles/claude-runtime/src/server.ts))
  on port 5355 under `tini`; the server starts as root and drops to uid 1001 before anything else runs.
- Runs the agent loop in [Hermes Agent](https://github.com/NousResearch/hermes-agent) (pinned commit,
  `/opt/hermes`): one `python -m grokky_agent` process per chat session
  ([`agent/grokky_agent/`](../dockerfiles/claude-runtime/agent/grokky_agent/)), which embeds Hermes' `AIAgent`
  and exchanges JSON lines with the server over stdin/stdout. See [Operations](#operations).
- The provider key reaches each agent process only in the first line of its stdin; Hermes calls the
  provider's model API directly. No key is in an agent's environment or in a file.
- Ships the `datagrok-*` skills ([`plugin/`](../dockerfiles/claude-runtime/plugin/)).
- Health check: `GET /health` → `{status: "ok"}`.

**`mcp-server`** ([`Dockerfile`](../dockerfiles/mcp-server/Dockerfile)):

- Runs `dist/index.js` (compiled from [`src/index.ts`](../dockerfiles/mcp-server/src/index.ts))
  on port 3003.
- Exposes Datagrok operations as MCP tools (HTTP/streamable transport).
- Auth is passed per-request via `x-user-api-key` and `x-datagrok-api-url` headers, which each session's
  agent sends on every call.

### Container Discovery (browser side)

[`ClaudeRuntimeClient.discover()`](../src/claude/runtime-client.ts) finds both containers at
runtime via `grok.dapi.docker`:

```typescript
const [runtimeContainers, mcpContainers] = await Promise.all([
  grok.dapi.docker.dockerContainers.filter('name = "grokky-claude-runtime"').list(),
  grok.dapi.docker.dockerContainers.filter('name = "grokky-mcp-server"').list(),
]);
this.ws = await grok.dapi.docker.dockerContainers.webSocketProxy(containerId, '/ws');
this.mcpServerUrl = `${grok.dapi.root}/docker/containers/proxy/${mcpContainerId}/mcp`;
```

The WebSocket and MCP HTTP traffic both flow through Datagrok's reverse proxy, so the browser never
talks to the containers directly.

Container code compiles **inside** the image, so `grok publish` / `npm run build` (browser bundle
only) do not update it — rebuild the image and redeploy after any change under `dockerfiles/`.

---

## Message Protocol

All communication between the browser and `claude-runtime` is JSON over a single WebSocket connection.
Message types are defined in [`dockerfiles/claude-runtime/src/types.ts`](../dockerfiles/claude-runtime/src/types.ts).

### Incoming messages (browser → runtime)

| `type` | Purpose | Key fields |
|--------|---------|------------|
| `user_message` | Send a prompt, start or resume a session | `sessionId`, `message`, `images?`, `apiKey`, `mcpServerUrl`, `outputSchema?`, `systemPromptMode?` (`datagrok`/`bash`/`none`), `model?`, `clientTools?`, `gates?`, `resumeExpected?`, `taskId?` |
| `abort` | Cancel an in-progress query | `sessionId` |
| `input_response` | Reply to an `input_request` (browser tool result or user answer) | `sessionId`, `requestId`, `value` |
| `sync_user_files` | Trigger a skills/knowledge sync | `apiKey`, `mcpServerUrl`, `scope?`, `packageName?` |
| `models` | Ask for the runtime's models: answered once one of its agents has started, or with the reason none could | — |

Notes on `user_message` fields:

- `images` — base64 attachments (`ImageAttachment[]`) prepended to the prompt content.
- `clientTools` — tool definitions declared by the browser for this turn (the static view-function
  meta-tools plus `datagrok_confirm`); their calls round-trip to the browser via
  `input_request`/`input_response`.
- `gates` — per-turn switches for the verify/grounding gates (both default on; used by benchmark arms).
- `taskId` — id of the queued admission task: the turn holds until the matching
  `Grokky:aiChatTurnTask` celery task claims it via `POST /task/claim` (see
  [`tasks.ts`](../dockerfiles/claude-runtime/src/tasks.ts), [`src/claude/queue-task.ts`](../src/claude/queue-task.ts)).

### Outgoing messages (runtime → browser)

| `type` | Purpose | Key fields |
|--------|---------|------------|
| `chunk` | Incremental text delta from Claude | `sessionId`, `content` |
| `tool_activity` | Human-readable summary of a tool call in progress | `sessionId`, `summary`, `name` (bare tool name) |
| `revision_start` | A gate (verifier/grounding) blocked the turn's Stop; a revision streams hidden — `final.revision` says whether it replaces the visible answer | `sessionId` |
| `final` | Completed response | `sessionId`, `content`, `structured_output?`, `unverified?`, `metrics?` (tokens/cost/turns/duration), `revision?` (`kept`/`replaced`) |
| `error` | Session error | `sessionId`, `message` |
| `session_reset` | The session record is gone but the browser sent `resumeExpected` — the turn was not run; the browser resends the prompt with the conversation transcript prepended | `sessionId` |
| `aborted` | Confirmed abort | `sessionId` |
| `queued` | Turn is waiting — behind the session's active query, or for its celery admission claim (re-sent every 60s while waiting) | `sessionId` |
| `input_request` | Runtime needs the browser to execute a tool or answer a question | `sessionId`, `requestId`, `toolName`, `input` |
| `sync_status` | Result of a `sync_user_files` request | `status`, `files?`, `message?` |
| `models` | The provider's model ids a turn may name, the one a turn without a model runs on, and the one for quick calls; or why no agent could start | `models`, `defaultModel`, `fastModel`, or `error` |

---

## Server-Side: claude-runtime

[`server.ts`](../dockerfiles/claude-runtime/src/server.ts) is a [Hono](https://hono.dev) server with a
WebSocket endpoint at `/ws`, `GET /health`, and `POST /task/claim` / `/task/release` for the celery
admission long-poll. [`session.ts`](../dockerfiles/claude-runtime/src/session.ts) runs the turns;
[`agent/pool.ts`](../dockerfiles/claude-runtime/src/agent/pool.ts) owns the agent processes and
[`agent/turn.ts`](../dockerfiles/claude-runtime/src/agent/turn.ts) drives one turn on one of them. Every
message shape — WebSocket protocol, turn request, agent commands and events — is in
[`types.ts`](../dockerfiles/claude-runtime/src/types.ts).

### Process model

```
tini (pid 1)
└── node dist/server.js                  drops from root to uid 1001 at start, so it is not dumpable
    ├── python -m grokky_agent           one per chat session: own process group, own HERMES_HOME
    │   │   stdin   init (credentials, settings), then turn / continue / cancel / browser_result
    │   │   stdout  events: ready, session, text, tool_start, tool_end, browser_call, stopped, …
    │   ├── Hermes AIAgent ──▶ the provider's model API (apiKey)
    │   ├── MCP client     ──▶ mcp-server (the user's Datagrok key in the headers)
    │   └── tool processes (terminal, …)
    └── python -m grokky_agent           warm spare, not yet bound to a session

/var/lib/grokky/hermes/<uuid>/           a live agent's HERMES_HOME; it stays, with its state.db, after the agent
                                         dies holding its session, and the next turn reloads it there
/var/lib/grokky/hermes/parked/<uuid>/    state.db of an evicted session, reloaded by its next turn
```

### Turn lifecycle

[`handleMessage()`](../dockerfiles/claude-runtime/src/session.ts) handles each `user_message`:

1. Chains the turn on the session's FIFO queue (`sessionChains`) — a second prompt for the same
   session waits behind the first, emitting `queued` heartbeats.
2. If the message carries a `taskId`, waits for the celery admission claim (the queue container's
   concurrency limit is the only cap on concurrent turns), failing open if no claim arrives.
3. Resolves the per-user directory (`/users/<id>`, the user verified with their key), which becomes
   the turn's working directory, and syncs the user's files (the first sync per user blocks; later
   turns serve what is on disk).
4. Picks the [system prompt](#system-prompt) for `systemPromptMode` and the turn's limits: 100 model
   calls and $2.50 in full mode, 10 calls and $0.25 in `bash`/`none`. Full mode also gets the gates
   below and `high` reasoning effort.
5. [`runTurn()`](../dockerfiles/claude-runtime/src/agent/turn.ts) answers `session_reset` if the browser
   expects to resume (`resumeExpected`) but the runtime holds neither a live process nor a parked session
   for it; the browser then replays its transcript once. Otherwise it takes the session's agent process
   from the pool (live, reloaded from its parked store, a warm spare, or a new one) and sends a `turn`
   command: text and images, system prompt, effort, the browser's model id, cwd, output schema, the browser's
   `clientTools`, the Datagrok MCP URL with the user's key in its headers, and the limits.
6. The agent ([`session.py`](../dockerfiles/claude-runtime/agent/grokky_agent/session.py)) prepares its
   `AIAgent` — the turn's model, or the provider's default (another model continues the same stored
   conversation; a model the provider does not offer ends the turn with an error, see
   [Credentials](#credentials)), this turn's
   view tools and MCP server (re-registered only when they change), thinking on or off — and runs the
   conversation. The hooks registered by the in-process plugin
   ([`hermes_plugin/`](../dockerfiles/claude-runtime/agent/grokky_agent/hermes_plugin/)) act for the turn:
   - `pre_tool_call` reports `tool_start`, applies the [access guard](#security-model), and makes a relative
     `write_file`/`patch` path absolute against the turn's directory;
   - `post_tool_call` reports `tool_end`, which feeds the gates;
   - `pre_api_request` clears the previous attempt's error, `post_api_request` counts tokens and cost,
     `api_request_error` keeps the provider's error text, and `pre_auxiliary_call` / `post_auxiliary_call`
     track and count history compression;
   - the `llm_request` middleware puts the turn's system prompt and output format on every model
     request, in the request's format. If it cannot, it sends a request the API refuses instead, so nothing
     of the conversation leaves the process.
7. `session.ts` forwards the agent's events (mapping below), racing each against a 90 s idle watchdog
   ([`watchdog.ts`](../dockerfiles/claude-runtime/src/watchdog.ts)). Browser round-trips in flight don't
   count, and the agent sends an `activity` heartbeat at most every 30 s while a tool or compression call
   runs or the model streams anything, thinking and tool-call arguments included (Hermes stamps every
   stream event on its activity clock). Hermes' other stamps (its waits, wait notices, reconnects) do not
   count, so a provider that sends nothing still trips the watchdog; during a history compression, only
   once Hermes gives up on the compression call (after 300 s at least, 600 s at most).
   A turn that goes silent gets an `error`, and its session is parked: the process tree is killed and the
   next turn reloads the history in a fresh process.

The agent's tools are the same in every mode: Hermes' native toolsets `file` (`read_file`, `write_file`,
`patch`, `search_files`), `terminal` (`terminal`, `process_manage`), `skills` (`skills_list`, `skill_view`;
`skill_manage` is blocked) and `web` (`web_search`, `web_extract`); the browser tools below; and the
Datagrok MCP server's tools, named `mcp__datagrok__*`, whose results (like the web tools') Hermes hands to
the model inside its `<untrusted_tool_result>` envelope.

The agent counts model calls and cost per turn, gate revisions and structured-output retries included; at a
limit the run stops, a call past it is refused before it is sent, and the turn ends in an `error` naming the
turn or cost limit.

### Browser round-trips: `input_request` / `input_response`

[`browser_tools.py`](../dockerfiles/claude-runtime/agent/grokky_agent/browser_tools.py) registers the
tools that run in the Datagrok tab as native Hermes tools, so their results reach the model as they
are. A call is checked against the tool's schema, then relayed: the agent emits `browser_call`, the
server sends an `input_request` and waits for the matching `input_response` (correlated by
`requestId`, so parallel calls don't cross wires), and the answer goes back to the agent as
`browser_result`. A call the browser does not answer in time returns an error result to the model.

- **`datagrok_exec`** (300s timeout) — JS executed in the tab against the live view; an optional
  inline `verify` assertion runs in the same round-trip.
- **`datagrok_verify`** (60s) — standalone re-read-the-state assertion.
- **`datagrok_show_entities`** (30s) — renders interactive entity cards in the chat.
- **View tools** (120s) — whatever the browser declared in `clientTools`: the static meta-tools
  `list_view_functions` / `list_view_widgets` / `get_view_function_result` / `call_view_function`
  (searching and invoking the current view's `getFunctions()` set) plus `datagrok_confirm`
  (see [`src/ai/view-tools.ts`](../src/ai/view-tools.ts)). A declared tool needs a plain name and a
  description; a name another tool already has is skipped.

`AskUserQuestion` round-trips the same way with a 60s answer timeout: an unanswered question returns
"did not answer in time" to the model instead of holding the turn open.

### Verification and grounding gates

Full-prompt turns (not `bash`/`none`) are gated in Node: each finished tool call arrives as the agent's
`tool_end` event (a call the guard refused is not reported), and the stop check runs when the agent
reports `stopped`. The browser can switch either gate off per turn via `gates`.

- **[`verify.ts`](../dockerfiles/claude-runtime/src/verify.ts)** enforces act → verify: any
  `datagrok_exec` (including failed ones), or a browser, view or MCP tool call that does not read, marks
  the turn as needing proof (reads are recognized fail-closed: the domain tools by their `op` argument,
  other tools by their name's verb). A passing verification clears it; a stop with proof pending is
  blocked up to 3 times, after which `final` carries `unverified: true` and the panel shows a
  "Not verified" warning.
- **[`grounding.ts`](../dockerfiles/claude-runtime/src/grounding.ts)** checks at stop whether an
  answer makes platform how-to claims without having opened a source this turn, and demands a
  grounded revision if so. Small talk and structured-output turns skip it.

When a gate blocks the stop, the runtime emits `revision_start` and sends the agent a `continue`
command carrying the gate's feedback. The revision streams *hidden*
([`stream-filter.ts`](../dockerfiles/claude-runtime/src/stream-filter.ts) suppresses its chunks);
the visible answer stays put until `final.revision` says whether the revision replaces it
(`replaced`) or the original stands (`kept`, signaled by the model answering `NO_REVISION`).

### Agent event → WebSocket message mapping

The agent writes one JSON event per stdout line; [`turn.ts`](../dockerfiles/claude-runtime/src/agent/turn.ts)
turns them into turn events, which `forwardEvent()` in [`session.ts`](../dockerfiles/claude-runtime/src/session.ts)
emits:

| Agent event                     | → WS message                                                        |
|---------------------------------|---------------------------------------------------------------------|
| `text`                          | `chunk` (through the fence/revision buffer)                         |
| `tool_start`                    | `tool_activity`                                                     |
| `tool_end`                      | — (fed to the gates)                                                |
| `browser_call`                  | `input_request`; the reply goes back to the agent as `browser_result` |
| `stopped`, gates pass           | `final` (with the turn's `metrics`)                                 |
| `stopped`, a gate blocks        | `revision_start`, and a `continue` command to the agent             |
| `limit`                         | `error` (turn/cost limit text)                                      |
| `error`                         | `error` (e.g. the provider's message)                               |
| `session_lost`                  | `session_reset` (the stored history could not be restored)          |
| `activity`, `cancelled`         | — (`activity` only keeps the watchdog fed)                          |

`ready`, `fatal` and `session` belong to the process's startup and bookkeeping (see
[Agent processes and sessions](#agent-processes-and-sessions)). In a structured-output turn the agent
reports text deltas as `activity`, and the answer arrives whole in `final`.

Tool results are not forwarded — the browser sees only activity summaries, generated by
[`toolSummary()`](../dockerfiles/claude-runtime/src/tool-summary.ts) from `toolFormatters` and, for the Datagrok
MCP tools, their `op`; view tools and unknown tools read as their name.

### System prompt

[`prompts.ts`](../dockerfiles/claude-runtime/src/prompts.ts) defines one prompt per
`systemPromptMode`, built once at startup so it stays byte-stable for prompt caching:

- `datagrok` (default): `DATAGROK_PROMPT`, the inlined `datagrok-exec` and `datagrok-entities` skill
  bodies, and an index of every skill in `plugin/skills/` (name and description) to load with
  `skill_view`.
- `bash`: a short prompt to run the given command with the terminal tool and output only its stdout. Its
  text does not stream: the answer is the final message, since models remark before the tool call anyway.
- `none`: no system prompt at all.

The agent's `llm_request` middleware
([`formats.py`](../dockerfiles/claude-runtime/agent/grokky_agent/formats.py), by the format Hermes sends the
request in) puts the turn's prompt where that format has it — Anthropic's `system` (one block with
`cache_control: ephemeral`), OpenAI Responses' `instructions`, or Chat Completions' leading system or
developer message — or removes it in `none` mode, so the mode can change from turn to turn within a
session. A Responses call keeps a one-line neutral instruction in `none` mode: Hermes sends none without
instructions and fills empty ones with its own identity. Full mode runs with thinking at `high` effort;
`bash` and `none` turns ask for no thinking and run without gates, which Hermes maps per model: off on Haiku
5.5, between tool calls only on Sonnet 5.5, Opus 5.5 and Fable 5.1 keep their adaptive thinking, which they
cannot turn off, and GPT gets reasoning effort `none`.

`DATAGROK_PROMPT` instructs the model to:

- Use `workspace/`-relative paths for the public repo clone (absolute `/workspace` is blocked by the access guard).
- Route capability-shaped requests through the `datagrok-*` skills (`skill_view`).
- Reach the current view's operations through the view meta-tools, preferring view
  functions over generic `datagrok_exec` code.
- Reach server-side entities through the per-domain MCP tools, discovering `op` parameters instead
  of guessing names.
- Ground answers in sources opened this turn (`workspace/help/INDEX.md` for product questions,
  skills for code), never in memory.
- Use `AskUserQuestion` only when the intent itself is unclear.

---

## Browser Side: ClaudeRuntimeClient

[`ClaudeRuntimeClient`](../src/claude/runtime-client.ts) is a singleton WebSocket client.

**API:**

| Method | Purpose |
|--------|---------|
| `send(sessionId, message, options?)` | Send a prompt; options: `outputSchema`, `systemPromptMode`, `model`, `taskId`, `clientTools` |
| `abort(sessionId)` | Cancel an in-progress query |
| `respondToInput(sessionId, requestId, value)` | Reply to an `input_request` |
| `query(message, options?)` | One-shot promise-based call with timeout; wraps send/subscribe/cleanup |
| `models()` | The runtime's model ids, default and fast model (asked once, then cached) |
| `syncUserFiles(scope?, packageName?)` | Trigger a skills/knowledge sync |
| `isResumable(sessionId)` | Whether the runtime holds this session (drives `resumeExpected`) |
| `expectResume(sessionId)` | Marks a session the runtime may still hold (a chat reopened from the history) |

**RxJS subjects** (one per outgoing message type): `onChunk`, `onToolActivity`, `onRevisionStart`,
`onFinal`, `onError`, `onSessionReset`, `onAborted`, `onInputRequest`, `onSyncStatus`, plus `onClose`
for the socket itself.

---

## Streaming to the Panel UI

The assistant UI is ONE singleton [`AIPanel`](../src/ai/panel.ts) (created by `initAIWindow()`),
used for every view. [`runPromptWithLifecycle()`](../src/ai/ui.ts) is the central routing function —
it intercepts `!`/`!!` command prefixes and usage limits, then `runClaudeStreaming()` in
[`src/ai/ui.ts`](../src/ai/ui.ts) connects `ClaudeRuntimeClient` events to the panel:

```
client.onChunk         → panel.updateStreaming(accumulated)           (live markdown update)
client.onToolActivity  → panel.updateStreaming(accumulated + status)  (tool badge)
client.onRevisionStart → delayed "Revising…" status
client.onFinal         → panel.finalizeStreaming(content)             (render; swap in revision if 'replaced')
client.onError         → grok.shell.error(msg) + cleanup
client.onAborted       → panel.clearStreaming() + cleanup
client.onSessionReset  → resend prompt with restored transcript prepended (once)
client.onInputRequest  → dispatch by toolName:
    datagrok_exec          → executeSingleBlock() (+ inline verify) → respondToInput(outcome)
    datagrok_verify        → runVerification() → respondToInput({passed, observed})
    datagrok_confirm       → approval card, held until a button is clicked
    datagrok_show_entities → renderEntityRefList() entity cards
    view tools             → viewFunctionTools() runner against the live grok.shell.v
    AskUserQuestion        → panel.showInputRequest() inline form
```

Every prompt gets a fresh workspace snapshot from
[`buildWorkspaceContext()`](../src/claude/exec-blocks.ts) (current view details and AI briefing,
widget briefings, all open views, all workspace tables) prepended browser-side, and full-mode turns
declare the static view-function meta-tools as `clientTools`. When the queued route is deployed,
`sendChatTurn()` ([`src/claude/queue-task.ts`](../src/claude/queue-task.ts)) first launches the
`aiChatTurnTask` celery task that holds the admission slot; the turn itself streams over the
WebSocket as usual.

Claude executes JS through the `datagrok_exec` tool *during* the turn (an `input_request`
round-trip via [`exec-blocks.ts`](../src/claude/exec-blocks.ts) `executeSingleBlock()`), so it
responds only after knowing the real outcome — responses are not post-processed for fenced code
blocks. If the executed code returns an `HTMLElement`, it is appended to the chat as a result
widget.

### Panel entry points

All are setup functions in [`src/ai/ui.ts`](../src/ai/ui.ts) called from `init()` in
[`src/package.ts`](../src/package.ts); the ribbon icons toggle the singleton panel rather than
creating per-view panels:

| Function | Context |
|----------|---------|
| `setupShellAIPanelUI()` | Shows the singleton panel in the AI window |
| `setupTableViewAIPanelUI()` | AI ribbon icon on table views |
| `setupScriptsAIPanelUI()` | AI ribbon icon on script views |
| `setupAIQueryEditorUI()` | Called by the core query editor — just reports whether AI is configured |
| `setupSearchUI()` | Global search bar — routes to `CombinedAISearchAssistant` |
| `setupAgentScriptsUI()` | Run button on file views under `MyFiles/agents/scripts/` |

---

## MCP Server Tools

[`dockerfiles/mcp-server/src/index.ts`](../dockerfiles/mcp-server/src/index.ts) registers tools using
`@modelcontextprotocol/sdk`. Tool calls are authenticated per-request via the headers each session's
agent sends (the user's key, which the runtime passes with every turn); all actual HTTP calls are made by
[`src/api-client.ts`](../dockerfiles/mcp-server/src/api-client.ts).

**One tool per domain, dispatching on `op`.** Each of the 34 operations used to be its own MCP
tool, so all 34 schemas sat in the model's prompt prefix on every turn — measured at 3,847
tokens/call. Collapsing them to five tools whose descriptions carry only the operation *names*
brought that to 1,061 tokens/call, and a call with no `op` returns the full parameter schema for
that domain, so nothing became undiscoverable. Registry:
[`src/ops.ts`](../dockerfiles/mcp-server/src/ops.ts).

| Tool | Operations |
|----------|-------|
| `datagrok_functions` | `list`, `get`, `call`, `list_scripts`, `list_queries`, `create_script`, `create_query` |
| `datagrok_files` | `list`, `download`, `upload` |
| `datagrok_projects` | `list`, `get`, `create`, `delete`, `search`, `list_recent`, `attach_entity`, `share`, `list_shares` |
| `datagrok_spaces` | `list`, `get`, `create`, `delete`, `create_subspace`, `list_children`, `add_entity`, `remove_entity`, `read_file`, `write_file`, `delete_file` |
| `datagrok_platform` | `whoami`, `list_users`, `list_groups`, `list_connections` |

Results are compact JSON, arrays are paged (50 by default, `offset`/`limit` to page), and
everything is capped — a tool result stays in the transcript and is re-read on every later hop, so
an unpaged list taxes the whole turn ([`src/format.ts`](../dockerfiles/mcp-server/src/format.ts)).

**Op names are load-bearing**: the runtime's verify gate decides whether a call mutated anything by
matching the *op* against read-verb prefixes (`list`/`get`/`search`/`read`/`download`/`find`/`whoami`).
An op named outside that convention silently changes whether verification is demanded, so
`mcp-server/test/ops.test.js` asserts every op starts with a known read or write verb.

---

## End-to-End Flow Example

Below is the full path for "Add a scatter plot of Height vs Weight":

```
1. User types in the AI panel → runPromptWithLifecycle() → UsageLimiter.tryCheckAndIncrement()
2. runClaudeStreaming() → client.ensureConnected(); buildWorkspaceContext() prepended
3. sendChatTurn() holds the celery admission slot (when deployed), then
   WS → {type:'user_message', sessionId, message, apiKey, mcpServerUrl, clientTools}

4. claude-runtime handleMessage():
   - agent/turn.ts runTurn(): {type:'turn', …} to the session's agent process; the agent's
     llm_request middleware puts DATAGROK_PROMPT on each model request

5. Claude reads JS API source (read_file on workspace/js-api/src/viewer.ts)
   → pre_tool_call (access guard) → agent event tool_start
   → {type:'tool_activity', summary:'Read workspace/js-api/src/viewer.ts'} → browser shows the badge

6. Claude calls datagrok_exec with the viewer code + inline verify assertion
   → agent event browser_call → {type:'input_request', toolName:'datagrok_exec', input:{code, verify}}
   → browser executeSingleBlock() adds the scatter plot to the live view,
     runVerification() re-reads view.viewers, respondToInput({success:true, verified:{passed:true}})
   → browser_result to the agent; post_tool_call → tool_end → the verify gate marks the action proven

7. Claude emits text deltas
   → agent event text → {type:'chunk', content:'Added a scatter plot…'}
   → browser onChunk → accumulated text, panel.updateStreaming()

8. The agent reports stopped, gates pass → {type:'final', content, metrics}
   → browser onFinal → panel.finalizeStreaming() renders the final markdown
```

---

## Key Files Reference

Runtime paths are relative to [`dockerfiles/claude-runtime/`](../dockerfiles/claude-runtime/).

| File | Purpose |
|------|---------|
| [`Dockerfile`](../dockerfiles/claude-runtime/Dockerfile) | Hermes stage (pinned), public repo in `/workspace`, `src/` via tsc, the agent package in `/app/agent` |
| [`entrypoint.sh`](../dockerfiles/claude-runtime/entrypoint.sh) | Execs the server under `env -i` with an allowlist |
| [`src/privileges.ts`](../dockerfiles/claude-runtime/src/privileges.ts) | The server's first import: drops from root to uid 1001 |
| [`src/server.ts`](../dockerfiles/claude-runtime/src/server.ts) | Hono server: `/ws` dispatch, `/health`, task claim/release, pool start, shutdown |
| [`src/session.ts`](../dockerfiles/claude-runtime/src/session.ts) | Per-session queue, file sync, gates, idle watchdog, browser round-trips, event forwarding, abort |
| [`src/agent/pool.ts`](../dockerfiles/claude-runtime/src/agent/pool.ts) | Agent processes: spawn with `init`, warm spares, LRU and idle eviction, parking and reload, process-tree kills |
| [`src/agent/turn.ts`](../dockerfiles/claude-runtime/src/agent/turn.ts) | One turn over the JSON-lines protocol: commands out, events in, gate continuation, cancel |
| [`src/types.ts`](../dockerfiles/claude-runtime/src/types.ts) | WebSocket protocol, turn request, `AgentCommand` / `AgentProcessEvent`, tool inputs |
| [`src/prompts.ts`](../dockerfiles/claude-runtime/src/prompts.ts) | System prompts, inlined skills, skills index (`buildSystemPrompt`) |
| [`src/tool-summary.ts`](../dockerfiles/claude-runtime/src/tool-summary.ts) | Tool-activity summaries |
| [`src/verify.ts`](../dockerfiles/claude-runtime/src/verify.ts) | Act → verify gate |
| [`src/grounding.ts`](../dockerfiles/claude-runtime/src/grounding.ts) | Grounding gate for platform how-to claims |
| [`src/stream-filter.ts`](../dockerfiles/claude-runtime/src/stream-filter.ts) | Fence-aware chunk filtering, hidden revision streaming |
| [`src/watchdog.ts`](../dockerfiles/claude-runtime/src/watchdog.ts) | Idle race, process-tree kills, orphan reaper |
| [`src/tasks.ts`](../dockerfiles/claude-runtime/src/tasks.ts) | Celery admission-task claim/release registry |
| [`src/constants.ts`](../dockerfiles/claude-runtime/src/constants.ts) | Paths (with their test overrides), the credential names handed to agents |
| [`agent/grokky_agent/__main__.py`](../dockerfiles/claude-runtime/agent/grokky_agent/__main__.py) | Agent entry: non-dumpable, `init`, plugin check, command loop |
| [`agent/grokky_agent/session.py`](../dockerfiles/claude-runtime/agent/grokky_agent/session.py) | The session on Hermes' `AIAgent`: hooks, limits, cancel, structured output, heartbeat, browser relay |
| [`agent/grokky_agent/hermes_plugin/`](../dockerfiles/claude-runtime/agent/grokky_agent/hermes_plugin/) | In-process Hermes plugin that registers the session's hooks and middleware |
| [`agent/grokky_agent/browser_tools.py`](../dockerfiles/claude-runtime/agent/grokky_agent/browser_tools.py) | Browser and view tools as native Hermes tools |
| [`agent/grokky_agent/guard.py`](../dockerfiles/claude-runtime/agent/grokky_agent/guard.py) | Access guard on tool arguments, relative edit paths |
| [`agent/grokky_agent/hermes.yaml`](../dockerfiles/claude-runtime/agent/grokky_agent/hermes.yaml), [`config.py`](../dockerfiles/claude-runtime/agent/grokky_agent/config.py) | Hermes' config for every session, written to `HERMES_HOME/config.yaml` with the skills directory, compression threshold and each model's context window and image input |
| [`agent/grokky_agent/models.yaml`](../dockerfiles/claude-runtime/agent/grokky_agent/models.yaml), [`models.py`](../dockerfiles/claude-runtime/agent/grokky_agent/models.py) | The model catalog: each model's family, context window, image input and prices (the cost limit) |
| [`agent/grokky_agent/formats.py`](../dockerfiles/claude-runtime/agent/grokky_agent/formats.py) | The request formats (Anthropic Messages, OpenAI Responses, Chat Completions): where the system prompt and output schema go, the refused request, the strict schema |
| [`agent/grokky_agent/hermes.py`](../dockerfiles/claude-runtime/agent/grokky_agent/hermes.py) | Every Hermes API the driver calls, typed |
| [`agent/grokky_agent/protocol.py`](../dockerfiles/claude-runtime/agent/grokky_agent/protocol.py) | The JSON-lines protocol (`AgentCommand` / `AgentProcessEvent` in `src/types.ts`) |
| [`agent/grokky_agent/structured.py`](../dockerfiles/claude-runtime/agent/grokky_agent/structured.py) | Output-schema check (`jsonschema`) |
| [`agent/grokky_agent/providers/`](../dockerfiles/claude-runtime/agent/grokky_agent/providers/) | `Provider` base (model credentials, format per family, test upstream, Hermes model metadata) and one module per provider (`AnthropicProvider`): credential checks, model ids, connection |
| [`dockerfiles/mcp-server/src/index.ts`](../dockerfiles/mcp-server/src/index.ts) | MCP tool registry |
| [`dockerfiles/mcp-server/src/api-client.ts`](../dockerfiles/mcp-server/src/api-client.ts) | HTTP client to Datagrok API |
| [`src/claude/runtime-client.ts`](../src/claude/runtime-client.ts) | Browser WebSocket client, RxJS subjects, `ClaudeRuntimeClient` singleton |
| [`src/claude/exec-blocks.ts`](../src/claude/exec-blocks.ts) | `executeSingleBlock()`, `runVerification()`, `buildWorkspaceContext()`, entity cards |
| [`src/claude/queue-task.ts`](../src/claude/queue-task.ts) | Queued-task admission (`aiChatTurnTask` holds the celery slot) |
| [`src/ai/view-tools.ts`](../src/ai/view-tools.ts) | View-function meta-tools declared as `clientTools` |
| [`src/ai/ui.ts`](../src/ai/ui.ts) | `runPromptWithLifecycle()`, streaming orchestration, panel setup functions |
| [`src/ai/panel.ts`](../src/ai/panel.ts) | `AIPanel` — the singleton chat panel (`StreamingPanel` interface) |

---

## Operations

How to configure, run, test, debug and roll back the `claude-runtime` container. Paths are relative to
[`dockerfiles/claude-runtime/`](../dockerfiles/claude-runtime/).

### Agent processes and sessions

See the [process model](#process-model).

- PID 1 is `tini` (reaps orphaned tool processes); it runs `entrypoint.sh`, which execs the server under
  `env -i` with an allowlist: the variables in [Configure](#configure), the credentials and a few fixed
  paths. The server's first import (`src/privileges.ts`) drops it from root to uid 1001.
- Each chat session gets one long-lived agent process (`src/agent/pool.ts`): `python -m grokky_agent` run
  with Hermes' venv Python, in a process group of its own, with cwd `/` (root-owned, so a stray relative
  path never lands in runtime state) and its own `HERMES_HOME` at `/var/lib/grokky/hermes/<uuid>`
  (generated `config.yaml`, a link to the plugin, Hermes' session store `state.db`, logs; kept off `/tmp`,
  which an agent may clear). Its environment holds only `PATH`, `HOME`/`HERMES_HOME`, `TERM`, `LANG`,
  `PYTHONUNBUFFERED` and `PYTHONPATH`.
- Startup: the first stdin line is `init` (the provider credentials and settings). The agent writes its
  config, loads the plugin, and answers `ready` with the provider's model ids, default and fast model, or
  `fatal` with the reason (a refused credential, the plugin not loading). An agent that is not ready within
  the ready timeout is disposed and its turn errors with "The assistant could not start: …".
- Warm spares are started ahead with Hermes loaded and an `AIAgent` built for the default model; a new chat
  takes a ready spare, or takes over one that is still starting rather than starting a second process.
- An idle session's process is evicted after the idle TTL (checked every minute), or earlier at the process
  cap when a new chat has no spare to take: spares go first, then the least recently used idle session. A busy session is never evicted, so the
  pool can run over the cap. Eviction kills the process tree and, once the process has exited, moves its
  `state.db` (with `-wal`/`-shm`) to `/var/lib/grokky/hermes/parked/<uuid>/`. The next turn starts a fresh
  process on that home, which first deletes everything in it but the session store, and continues the
  stored session. Past the session cap (live plus parked) the oldest parked session is deleted. A session
  whose stored history is gone answers `session_reset`; a reload that fails to start keeps the session
  parked, so the turn errors and the next one retries. A process that dies holding a session (killed, out
  of memory, crashed) keeps its session store where it is, and the next turn reloads it; a reloaded process
  holds the stored session until its agent has opened it, so a turn that fails before that (an unusable
  output schema) or is stopped while the process starts loses nothing either.
- The browser saves each chat's session id with the conversation: a chat reopened from the history, after a
  page reload too, resumes that session while the runtime still has it, and otherwise gets
  `session_reset` and sends the turn again with the saved transcript.
- Abort: the server sends `cancel`; the agent interrupts Hermes and fails its pending browser calls with
  `aborted`. Once the turn has settled, the server kills everything under the agent process, since Hermes
  does not stop what a tool started in a session of its own. The next turn continues the conversation,
  which keeps the interrupted request.
- Disposing a process kills its whole descendant tree (read from `/proc`), then its process group. SIGTERM
  or SIGINT disposes every agent and exits; startup deletes whatever a previous run left under
  `/var/lib/grokky/hermes`, parked sessions included. An orphan reaper (every 30 s) kills any
  `-m grokky_agent` process the pool does not track, with its tree.
- Hermes compresses the history once the estimated prompt passes the compression threshold, keeping the
  recent messages verbatim; the compression call counts toward the turn's cost.
- The Hermes config (`agent/grokky_agent/hermes.yaml`) turns Hermes' own additions off — nudges and
  guidance, tool search, connectors, loop guardrails, approvals, title generation, the remote model
  catalog, borrowed logins, OpenRouter's response cache — because the runtime owns the prompt and the loop's
  policy. `AIAgent` runs without context files (no `AGENTS.md`/`CLAUDE.md` from the workspace), memory or
  background review.
  `agent.api_max_retries: 4` is the number of attempts per model call.
- Web tools: with no web backend configured, Hermes' `web_search` and `web_extract` rotate through the
  keyless free tiers of Exa, Parallel, Firecrawl and Keenable (`web.keyless_fallback`, on by default), so
  queries and URLs go to those vendors.
- The Datagrok MCP server is connected when a turn first carries it. The connect, the MCP handshake and
  the tool listing are bounded by 30 s (Hermes' default is 60 s), so a slow server cannot hold the turn
  past the idle watchdog. A server that could not be reached is dropped, and the next turn connects again:
  while it refuses connections that costs each turn Hermes' three quick retries (about 7 s); one that does
  not answer costs the 30 s, then Hermes waits 30 s (doubling) before it tries again.

### Configure

| Env var (server) | Default | Meaning |
|---|---|---|
| `HERMES_MAX_PROCESSES` | 12 | live agent processes, spares included; at the cap the least recently used idle session is evicted (logged `agent: evicting session … (LRU, cap 12)`) |
| `HERMES_IDLE_TTL_MS` | 900000 | an idle session's process is evicted after this (logged `… (idle N min)`); its next turn reloads it |
| `HERMES_MAX_SESSIONS` | 200 | live plus parked sessions; past it the oldest parked session is deleted (logged `agent: dropped the stored state of session …`) and gets `session_reset` |
| `HERMES_WARM_SPARES` | 2 | agent processes started ahead of need |
| `HERMES_READY_TIMEOUT_MS` | 30000 | an agent that has not answered `ready` by then is disposed; its turn errors |
| `HERMES_COMPRESSION_THRESHOLD_TOKENS` | 160000 | history compression trigger (estimated prompt tokens) |
| `HERMES_LOG_STDERR` | unset | any value forwards every agent stderr line to the runtime log |
| `TEST_UPSTREAM` | unset | test stubs and loggers only: replaces the scheme, host and port of every model URL of the provider (the path stays). It may name http(s) to `127.0.0.1`, `localhost`, `[::1]` or `host.docker.internal`, or https to the provider's own host on another port; with any other value every agent refuses to start |
| `DATAGROK_API_URL` | injected by Datlas | the URL user identity is verified against; without it, the first client-supplied URL is trusted (dev only) |

`entrypoint.sh` passes each of these variables through `env -i`, so
`docker run -e HERMES_LOG_STDERR=1 …` works. The agent homes
(`/var/lib/grokky/hermes`), skills (`/app/plugin`), agent package (`/app/agent`) and agent Python
(`/opt/hermes/.venv/bin/python`) are fixed by the image; a server run outside it overrides them with
`GROKKY_HERMES_HOMES`, `GROKKY_PLUGIN_DIR`, `GROKKY_AGENT_DIR` and `GROKKY_AGENT_PYTHON`
(`src/constants.ts`). `HERMES_RUNTIME_DIR` is Hermes' own package-manager store, not a runtime setting.

Hermes' `config.yaml` is written by each agent at startup from `agent/grokky_agent/hermes.yaml` (with the
skills directory, the compression threshold, and each of the provider's models' context window and image
input as `model_overrides`); change settings there. Memory: an idle agent holds about 250 MB RSS, so the 12
processes, spares included, fit the 8 GB limit in `container.json`.

### Credentials

The provider settings and key are the `claude-runtime` container's credentials: the keys and values of the
container's Credentials dialog, which the platform passes to the container as environment variables. They
are not Grokky's package credentials, which reach only the browser; the same names are to move there once
the platform passes package credentials to a container. The container's identity includes the package
version, so they probably have to be entered again after a Grokky release (inferred from the platform's
code, not tested).

| Credential | Meaning |
|---|---|
| `provider` | The model provider: absent or `Anthropic`, or another provider below |
| `apiKey` | The provider's key; required |
| `baseUrl` | The provider's address, where the admin chooses it (a resource, a workspace, an endpoint) |
| `region` | The provider's region, where its hosts are regional |
| `apiFormat` | A custom endpoint's request format: `anthropic-messages`, `openai-responses` or `openai-chat-completions` |
| `model` | The model a turn runs on when the browser names none: one of the provider's model ids, or `<id>=<catalog model>` for a model the admin names, such as a deployment |
| `fastModel` | The model for quick calls (search ranking, query matching, the benchmark's judge), in the same forms; by default the provider's fast model, else `model` |

The server takes these (and `TEST_UPSTREAM`, see [Configure](#configure)) out of its environment at
startup and sends them to each agent in its `init` line (`PROVIDER_CREDENTIALS` in `src/constants.ts`); the
agent's provider (`agent/grokky_agent/providers/`) checks them and passes the key to Hermes in memory. An
agent refuses to start, so every turn ends in an `error` "The assistant could not start: …", when its
provider refuses the credentials (the reasons are in each provider's subsection), when `provider` names no
provider of the runtime, or when `TEST_UPSTREAM` is not an allowed host. The server itself stays up; the
log shows the reason as `agent: warm spare failed to start: …`. A key the provider rejects ends the turn in
an `error` with the provider's message; a 401 adds a hint naming the credential to check. An answer the
provider's safety filter stops mid-stream ends in an `error` that says so, not with the cut-off text.

**Models and prices.** The runtime ships a model catalog,
[`models.yaml`](../dockerfiles/claude-runtime/agent/grokky_agent/models.yaml): each catalog model's family
(`claude` or `gpt`), context window, image input, and prices in USD per million tokens (input, output, cache
reads and writes, and the long-context rates where a long prompt is billed at other rates), with the
pricing pages they come from. Each provider maps its model ids to catalog models (`model_ids` of its class)
and names its default and fast model; `model` and `fastModel` change those two, and map a model the admin
names to the catalog model it runs. Claude is called in Anthropic's Messages format and GPT in OpenAI's
Responses format, unless a provider says otherwise. A turn runs only on one of the provider's models: the
browser names one or none, and any other model ends the turn with an error, so the cost limit always has a
price. The agent prices every model call, compression included, at its catalog model's rates; Hermes gets
each model's context window and image input, and a turn with images on a model without image input ends
with an error. The agent's `ready` event reports the provider's model ids with its default and fast model,
and the browser asks for them (`models`, `ClaudeRuntimeClient.models()`): Help and Execute run on the
default, search ranking and query matching on the fast model, and the benchmark offers the list.

Each provider's specifics follow in a subsection of its own. A provider is a `Provider` subclass
(`providers/base.py`) in a module of `providers/`, which registers it under its `name`.

#### Anthropic

`provider` absent or `Anthropic`; `apiKey` an Anthropic Console API key (`sk-ant-api…`). Hermes calls the
Messages API at `api.anthropic.com` with the key in `x-api-key`. Models: `claude-haiku-5-5` (fast),
`claude-sonnet-5-5` (default), `claude-opus-5-5` and `claude-fable-5-1`. An agent refuses a missing key, a
key with whitespace or control characters, and a Claude subscription token (an `sk-ant-` key other than
`api`, `usr`, `svc` or `admin`, a JWT or `cc-…`), which Hermes would send to `api.anthropic.com` as a
Claude Code login. Claude subscription logins are not supported: Anthropic allows them only in its own
applications, such as Claude Code, never in a third-party product
([Claude Code legal and compliance](https://code.claude.com/docs/en/legal-and-compliance#authentication-and-credential-use)).
The check applies to whatever key goes to `api.anthropic.com` (`Provider.refuse_subscription_token`), and
only there: other hosts take JWTs (Entra ID, Databricks OAuth).

#### Azure

`provider` `Azure`; `apiKey` a key of the resource (in the Foundry portal, a deployment's Details tab);
`baseUrl` the resource without a path: `https://<resource>.services.ai.azure.com` for a Microsoft Foundry
resource, which serves Claude and GPT, or `https://<resource>.openai.azure.com` for an Azure OpenAI resource,
GPT only; `model` (required) and `fastModel` the deployments to use, each as `<deployment>=<catalog model>`,
or as the deployment name alone when that is a catalog model id, as Azure names a deployment by default
(`claude-sonnet-5-5`, `gpt-6-sol`). The provider offers only these deployments; without `fastModel` the
default model is also the fast one. A Claude deployment's name must contain its model id and no dot
(`claude-opus-5-5-eu=claude-opus-5-5`): Hermes picks a Claude request's thinking and sampling settings by the
model name, so under another name Opus 5.5 and Fable 5.1 would get thinking with a token budget, which they
refuse, and it turns the dots of a `claude-…` name into hyphens. Of the catalog's models Azure deploys Claude
Opus, Sonnet and Haiku 5.5 and, in preview on Anthropic's infrastructure, Fable 5.1; GPT 6.1 Sol, 6 Astra,
Sol and Luna, and 5.6 Sol, Terra and Luna
([Claude in Foundry](https://learn.microsoft.com/en-us/azure/foundry/foundry-models/how-to/use-foundry-models-claude),
[models sold by Azure](https://learn.microsoft.com/en-us/azure/foundry/foundry-models/concepts/models-sold-directly-by-azure)).

Claude is called through the Messages API at `<baseUrl>/anthropic/v1/messages` with the key in `x-api-key`;
GPT through the Responses API of the
[v1 API](https://learn.microsoft.com/en-us/azure/foundry/openai/api-version-lifecycle) at
`<baseUrl>/openai/v1/responses`, without `api-version`, with the key as Bearer: OpenAI's client sends it so,
as in Microsoft's v1 examples with a resource key. Hermes would send a key to any `azure.com` Anthropic
address as Bearer, which Azure documents for Entra ID tokens only; `hermes.new_agent` turns Hermes' list of
such hosts off, so a string key always goes in `x-api-key` and only a callable one as Bearer. To Azure's
hosts Hermes also adds `api-version=2025-04-15` and the 1M-context beta to Claude calls, and replays less
encrypted reasoning on GPT calls; a stub sees this only when it answers as the resource
(`TEST_UPSTREAM=https://<resource>.services.ai.azure.com:<port>`).

An agent refuses a missing key or one with whitespace or control characters, a `baseUrl` of any other form
(a path, a port, another host or cloud, http), a Claude deployment on an `openai.azure.com` resource or under
another name, a missing `model`, and a deployment mapped to a model outside the catalog. A 401 adds a hint
to check `apiKey`; a 404, Azure's answer to an unknown deployment, one to check `baseUrl`, `model` and
`fastModel`. Microsoft Entra ID sign-in is not supported yet, so a resource with key access turned off
refuses every call.

### Build and run

```bash
cd packages/Grokky/dockerfiles/claude-runtime
docker build -t grokky-claude-runtime:<version> .
docker run --rm -p 5355:5355 -e apiKey grokky-claude-runtime:<version>   # apiKey exported in this shell
curl -s localhost:5355/health          # {"status":"ok"}
```

Hermes is pinned by `HERMES_COMMIT` in the Dockerfile's `hermes` stage and built with Hermes' own
package manager (CPython and uv sha-verified from `pm/lock.json`) with the `anthropic` and `mcp` extras.
The commit is recorded in the image label `ai.datagrok.hermes.commit` and in `/opt/hermes/COMMIT`. To
move the pin: change `HERMES_COMMIT`, rebuild, and run the BDD suite (see [Test](#test)) on a stand with
the new image.

### Security model

- The server and every agent run as uid 1001. The image refuses to start as any user but root
  (`entrypoint.sh`), because a server started as uid 1001 could not make itself non-dumpable.
- The server holds the provider credentials and every connected user's Datagrok key in its memory and
  sockets. It starts as root and drops to uid 1001 (gid 2001, no supplementary groups) before anything
  else runs (`src/privileges.ts`); a process whose uid changed after exec is not dumpable, so the agents'
  tools cannot read the server's `/proc/<pid>/{environ,mem,fd}`. It takes the provider credentials out of
  its environment at startup (`PROVIDER_ENV`, `src/constants.ts`): a child it runs (git, for the workspace
  sync) is dumpable again after its exec, and would otherwise inherit them.
- An agent gets the provider credentials only in its `init` line, never in its environment or a file,
  and keeps the key in memory. Before it imports the rest of the driver or reads stdin, it makes itself
  non-dumpable (`prctl(PR_SET_DUMPABLE, 0)` in `__main__.py`), so other agents cannot read its memory,
  environment or fds. The
  agent's stdin/stdout/stderr are socketpairs (Node `pipe` stdio), which cannot be reopened through
  `/proc/<pid>/fd`.
- The user's Datagrok key reaches the agent in each `turn` command as the MCP headers (`x-user-api-key`,
  and `Authorization` for the datlas container proxy).
- The provider key and the conversation go only to the provider's own hosts: a fixed address, or the
  admin's `baseUrl` under the provider's host rule. `TEST_UPSTREAM` can send them only to a local host, or to
  the provider's own host on another port. Nothing else in the image filters the agents' traffic.
- Hermes never borrows a Claude Code login (`auth.adopt_external_logins: false`), and logs at `WARNING`:
  `INFO` would write the start of every prompt into `agent.log`, which other agents can read.
- The access guard (`agent/grokky_agent/guard.py`) refuses tool calls that name another user's
  `/users/<id>`, a `..` segment, an absolute `/workspace` (a user's session works in its own
  `workspace/`), the agent homes (except the session's own `cache/spillover/`, where Hermes leaves a large
  tool result for `read_file`), a `state.db`, a `~` path in a file tool, or `/proc` beyond `/proc/self`
  and a few system files, and it refuses `skill_manage`. It is a pattern check on tool arguments: defence
  in depth, not isolation.
- Residual risks, all from the shared uid:
  - Agents of different users can read and write each other's files: user directories, live agent homes
    and parked session stores. A reload keeps only `state.db*` from the stored home and the agent writes
    its config fresh, so a planted `config.yaml`, `.env` or login never reaches a reloaded agent; a write
    into a live agent's home, or into a stored `state.db`, is not prevented.
  - The provider key is protected by the non-dumpable flag, not by a uid boundary.
  - An agent is dumpable between its exec and the `prctl` call. It has not read the credentials yet, but a
    same-uid process that attaches with `ptrace` in that window (the container has no Yama restriction)
    stays attached and can read them later.
  - `/workspace`, which every session reads (the help index the grounding gate sends the model to
    included), is writable by every agent.
  - A tool process that detaches into a new session of its own and its parent exits is no longer under the
    agent, so the cleanup after a stopped turn does not reach it; it lives until the container stops.
  - Any agent can signal the server and the other agents (same uid), so a hostile one can end every session
    in the container until it restarts.

### Test

```bash
npm ci
npm test              # tsc, then the unit tests
npm run lint
npm run lint:agent    # ruff, then mypy (strict; pyproject.toml), on the agent driver
```

`lint:agent` needs `python3` 3.11+ with the image's `jsonschema` and `ruamel.yaml`, plus `ruff`, `mypy` and
`types-jsonschema`. The AI panel's acceptance suite is the BDD feature set in `packages/Grokky/bdd`
(`GROKKY_BDD_LIVE=1 DATAGROK_SERVER=<alias> npx grok-bdd run --workers 1`), run against a stand with this
package and its runtime image deployed.

### Debug

- Runtime log lines: `agent[<home8>]: dispose (<reason>)`, `agent: evicting session <sid> (LRU|idle …|turn …
  went silent)`, `agent: reloaded session <sid> in N ms`, `agent: dropped the stored state of session …`,
  `agent[…]: exited unexpectedly (<code>); last stderr:`, `agent: warm spare failed to start: <reason>`,
  `agent[…]: a stdout line is not an event` (something printed on the event channel),
  `watchdog[<sid>]: no agent event for 90s — killing turn`, `reaper: SIGKILL orphaned …`,
  `task[<sid>]: no claim for … — failing open`, and per turn `verify[<sid>]: {…}` / `grounding[<sid>]: …`.
- Agent stderr: the last 40 lines of each process are kept for error messages; `HERMES_LOG_STDERR=1`
  forwards every line as `agent[<home8>] …`. The configured `apiKey` and key-shaped text (`sk-…`) are masked
  there and in the agents' error messages. Hermes' own logs
  (`logs/agent.log`, `logs/errors.log`, `WARNING` and up) are in the agent's home while the process lives;
  the home is deleted when the process ends.
- Processes in the container (the image has no `ps`):
  `docker exec <c> sh -c 'for p in /proc/[0-9]*; do printf "%s %s\n" ${p#/proc/} "$(tr "\0" " " <$p/cmdline)"; done'`.
  Expect `tini` (pid 1, root), then `node /app/dist/server.js` and one `python -m grokky_agent` per live
  session plus the spares, all as uid 1001. A non-dumpable process's `environ`, `mem` and `fd` need
  `CAP_SYS_PTRACE`, even for root: `docker exec --privileged -u 0`.
- What the model receives: the middleware sets the mode's prompt and, for `outputSchema` turns, the output
  format, in the request's format (`formats.py`): Anthropic's `system` and `output_config.format`; OpenAI
  Responses' `instructions`, `text.format` and `store: false` (every call resends the history, encrypted
  reasoning included); Chat Completions' leading system message and `response_format`. The model id and
  thinking come from the agent's `AIAgent`.
- A model request with `model: grokky-refused` means the middleware could not vouch for it (no turn in
  progress, or the patch failed); the API rejects it and nothing of the conversation went out.
- Structured output errors: "The output schema cannot be used: …" (before any model call) means the
  schema is not a JSON object or not a valid schema of its draft (draft-07 when it declares none).
  Anthropic's format gets the schema reduced to its strict grammar, OpenAI's formats get it as sent and not
  strict (strict mode needs every property required); the reply is validated against the original with
  `jsonschema`, and a mismatch is sent back to the model, up to 5 attempts: "… after 5 attempts" means
  every reply failed. The check runs in the session's own process, so a backtracking `pattern` stalls only
  that session, which the idle watchdog ends.
- Known limitation: Hermes' `search_files` with `target: "files"` matches file names only (a directory in
  the pattern is dropped) and also returns every hidden directory under `path` (`.git`, `.vscode`, …). The
  full-mode prompt tells the model; `none` mode shows it as is.

### Roll back

Redeploy the last SDK runtime image (the registry tag live before the migration) by pointing
`container.json` or the package version back at it and republishing. Roll back the Grokky package with
it: the SDK runtime reported a missing or rejected key with `auth_required`, which the current browser
package no longer handles, so an SDK image behind the current package leaves the panel waiting. No data
migration is involved; neither runtime keeps sessions across a container restart, so a rollback costs one
`session_reset` per open chat.
