# Claude Code Flow

This document describes how Claude Code is deployed and integrated into Datagrok: how the server-side
containers are set up, how messages travel between the browser and the containers, and how Claude's
responses are turned into live Datagrok UI.

---

## Architecture Overview

Three components work together:

```
Browser (TypeScript)
  └── ClaudeRuntimeClient  ──WebSocket──▶  claude-runtime container  ──MCP──▶  mcp-server container
                                                  │                                    │
                                          @anthropic-ai/claude-agent-sdk          grok.dapi.*
                                          (streams SDK events)                (Datagrok API over HTTP)
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

- Base image: `node:22-slim`.
- Clones the public Datagrok repo into `/workspace` at build time — this gives Claude read access to
  the full JS API source, packages, samples, and documentation without needing to fetch them at runtime.
- Runs `dist/server.js` (compiled from [`src/server.ts`](../dockerfiles/claude-runtime/src/server.ts))
  on port 5355.
- Model provider is configured via injected package credentials, translated to SDK env vars by
  `applyProviderConfig()` in [`server.ts`](../dockerfiles/claude-runtime/src/server.ts): Anthropic
  (`apiKey`), Amazon Bedrock (`awsBearerToken` or key pair), Microsoft Foundry
  (`foundryResource` + `foundryApiKey`), or a Claude subscription.
- Ships a local Claude Code plugin ([`plugin/`](../dockerfiles/claude-runtime/plugin/)) with the
  `datagrok-*` skills; attached for full-prompt sessions only.
- Health check: `GET /health` → `{status: "ok"}`.

**`mcp-server`** ([`Dockerfile`](../dockerfiles/mcp-server/Dockerfile)):

- Runs `dist/index.js` (compiled from [`src/index.ts`](../dockerfiles/mcp-server/src/index.ts))
  on port 3003.
- Exposes Datagrok operations as MCP tools (HTTP/streamable transport).
- Auth is passed per-request via `x-user-api-key` and `x-datagrok-api-url` headers forwarded from
  the claude-runtime on every call.

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
only) do not update it — rebuild the image and redeploy after any change under `dockerfiles/`. The
subscription re-auth flow (`auth_*` messages below) is walked through in the package
[`CLAUDE.md`](../CLAUDE.md).

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
| `auth_start` | Begin subscription re-auth (spawns `claude auth login`) | — |
| `auth_code` | Send the pasted OAuth code to the CLI's stdin | `code` |

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
| `tool_activity` | Human-readable summary of a tool call in progress | `sessionId`, `summary`, `name?` (bare tool name; absent for progress-only activity) |
| `revision_start` | A gate (verifier/grounding) blocked the turn's Stop; a revision streams hidden — `final.revision` says whether it replaces the visible answer | `sessionId` |
| `final` | Completed response | `sessionId`, `content`, `structured_output?`, `unverified?`, `metrics?` (tokens/cost/turns/duration), `revision?` (`kept`/`replaced`) |
| `error` | Session error | `sessionId`, `message` |
| `session_reset` | The session record is gone but the browser sent `resumeExpected` — the turn was not run; the browser resends the prompt with the conversation transcript prepended | `sessionId` |
| `auth_required` | The provider rejected the credentials mid-turn; the panel shows the re-auth widget | `sessionId` |
| `aborted` | Confirmed abort | `sessionId` |
| `queued` | Turn is waiting — behind the session's active query, or for its celery admission claim (re-sent every 60s while waiting) | `sessionId` |
| `input_request` | Runtime needs the browser to execute a tool or answer a question | `sessionId`, `requestId`, `toolName`, `input` |
| `sync_status` | Result of a `sync_user_files` request | `status`, `files?`, `message?` |
| `auth_url` | OAuth URL for the browser to open during re-auth | `url` |
| `auth_done` | Subscription re-auth succeeded | — |
| `auth_error` | Re-auth failed or none in progress | `message` |

---

## Server-Side: claude-runtime

[`server.ts`](../dockerfiles/claude-runtime/src/server.ts) is a [Hono](https://hono.dev) server with a
WebSocket endpoint at `/ws` plus `POST /task/claim` / `/task/release` for the celery admission
long-poll.

### Turn lifecycle

[`handleMessage()`](../dockerfiles/claude-runtime/src/session.ts) handles each `user_message`:

1. Chains the turn on the session's FIFO queue (`sessionChains`) — a second prompt for the same
   session waits behind the first, emitting `queued` heartbeats.
2. If the message carries a `taskId`, waits for the celery admission claim (the queue container's
   concurrency limit is the only cap on concurrent turns), failing open if no claim arrives.
3. Looks up the session record (in-memory LRU map, max 200 entries) and passes its SDK id as
   `resume`. After an abort, the next turn *forks* off the last clean assistant uuid instead,
   dropping the aborted turn from history.
4. Fires a background user-file sync and resolves the per-user directory (`/users/<id>`), which
   becomes the SDK `cwd`.
5. Calls `query()` from `@anthropic-ai/claude-agent-sdk` with the options built by
   [`buildOptions()`](../dockerfiles/claude-runtime/src/query-options.ts):
   - The [system prompt](#system-prompt) for the message's `systemPromptMode`.
   - Built-in tools: `Read`, `Glob`, `Grep`, `Edit`, `Write`, `Bash`, `WebSearch`, `WebFetch`,
     `AskUserQuestion`, `Skill` (`permissionMode: 'acceptEdits'`).
   - `mcpServers`: `datagrok` (HTTP, the mcp-server container), `datagrok-browser` (in-process,
     tools execute in the tab), and `datagrok-view` when the browser declared `clientTools`.
   - `model: model ?? sonnet`; full-prompt turns get thinking and high effort, `bash`/`none` are
     minimal. Per-mode `maxTurns`/`maxBudgetUsd` limits.
   - A `PreToolUse` access-guard hook confining file access to the user's own directory, and the
     gate hooks below.
6. Iterates the SDK event stream, forwarding each event via `forwardEvent()`. A watchdog bounds
   silence between events (90s; browser round-trips in flight don't count) — a wedged CLI is
   killed rather than hanging the session ([`watchdog.ts`](../dockerfiles/claude-runtime/src/watchdog.ts)).

### Browser round-trips: `input_request` / `input_response`

Tools that must run in the Datagrok tab are declared to the SDK as in-process MCP servers
([`query-options.ts`](../dockerfiles/claude-runtime/src/query-options.ts)); their handlers emit an
`input_request` and suspend until the matching `input_response` arrives (correlated by `requestId`,
so parallel calls don't cross wires):

- **`datagrok_exec`** (300s timeout) — JS executed in the tab against the live view; an optional
  inline `verify` assertion runs in the same round-trip.
- **`datagrok_verify`** (60s) — standalone re-read-the-state assertion.
- **`datagrok_show_entities`** (30s) — renders interactive entity cards in the chat.
- **View tools** (120s) — whatever the browser declared in `clientTools`: the static meta-tools
  `list_view_functions` / `list_view_widgets` / `get_view_function_result` / `call_view_function`
  (searching and invoking the current view's `getFunctions()` set) plus `datagrok_confirm`
  (see [`src/ai/view-tools.ts`](../src/ai/view-tools.ts)).

`AskUserQuestion` is intercepted in `canUseTool` and round-tripped the same way, but with a 60s
answer timeout: an unanswered question denies the tool call instead of holding the turn open.

### Verification and grounding gates

Full-prompt turns (not `bash`/`none`) are gated by a `PostToolUse`/`Stop` hook pair; the browser can
switch either off per turn via `gates`.

- **[`verify.ts`](../dockerfiles/claude-runtime/src/verify.ts)** enforces act → verify: any
  `datagrok_exec` (including failed ones) or non-read MCP call marks the turn as needing proof
  (reads are recognized fail-closed by the domain tools' `op` argument). A passing verification
  clears it; `Stop` blocks the turn from ending while proof is pending, up to 3 times, after which
  `final` carries `unverified: true` and the panel shows a "Not verified" warning.
- **[`grounding.ts`](../dockerfiles/claude-runtime/src/grounding.ts)** checks at Stop whether an
  answer makes platform how-to claims without having opened a source this turn, and demands a
  grounded revision if so. Small talk and structured-output turns skip it.

When a gate blocks Stop, the runtime emits `revision_start` and streams the revision *hidden*
([`stream-filter.ts`](../dockerfiles/claude-runtime/src/stream-filter.ts) suppresses its chunks);
the visible answer stays put until `final.revision` says whether the revision replaces it
(`replaced`) or the original stands (`kept`, signaled by the model answering `NO_REVISION`).

### SDK event → WebSocket message mapping

[`forwardEvent()`](../dockerfiles/claude-runtime/src/session.ts) maps SDK events:

| SDK event type               | → WS message                                     |
|------------------------------|--------------------------------------------------|
| `assistant` (tool_use block) | `tool_activity`                                  |
| `assistant` (auth failure)   | flags the session → `auth_required` at result    |
| `stream_event` (text_delta)  | `chunk` (filtered through the fence/revision buffer) |
| `tool_progress`              | `tool_activity`                                  |
| `tool_use_summary`           | `tool_activity`                                  |
| `result` (success)           | `final` (with `metrics`, commits the resume point) |
| `result` (error)             | `error`                                          |

Tool results are not forwarded — the browser sees only activity summaries, generated by
[`toolSummary()`](../dockerfiles/claude-runtime/src/query-options.ts) from the formatters in
`toolFormatters` and `mcpFormatters`.

### System prompt

[`prompts.ts`](../dockerfiles/claude-runtime/src/prompts.ts) defines one prompt per
`systemPromptMode`: `bash` (a minimal "output only stdout" prompt), `none` (empty — pure Claude
Code), and the full `DATAGROK_PROMPT` (default), which instructs Claude to:

- Use `workspace/`-relative paths for the public repo clone (absolute `/workspace` is hook-blocked).
- Route capability-shaped requests through the `datagrok-*` plugin skills (Skill tool).
- Reach the current view's operations through the `datagrok-view` meta-tools, preferring view
  functions over generic `datagrok_exec` code.
- Reach server-side entities through the per-domain MCP tools, discovering `op` parameters instead
  of guessing names.
- Ground answers in sources opened this turn (`workspace/help/INDEX.md` for product questions,
  skills for code), never in memory.
- Use `AskUserQuestion` only when the intent itself is unclear.

The `datagrok-exec` and `datagrok-entities` skill bodies are inlined into the prompt at startup;
the prefix is built once so it stays byte-stable for prompt caching.

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
| `syncUserFiles(scope?, packageName?)` | Trigger a skills/knowledge sync |
| `startAuth()` / `sendAuthCode(code)` | Drive the subscription re-auth flow |
| `isResumable(sessionId)` | Whether the runtime holds this session (drives `resumeExpected`) |

**RxJS subjects** (one per outgoing message type): `onChunk`, `onToolActivity`, `onRevisionStart`,
`onFinal`, `onError`, `onSessionReset`, `onAborted`, `onInputRequest`, `onSyncStatus`, `onAuthUrl`,
`onAuthDone`, `onAuthError`, `onAuthRequired`, plus `onClose` for the socket itself.

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
client.onAuthRequired  → auth renewal widget
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
`@modelcontextprotocol/sdk`. Tool calls are authenticated per-request via headers forwarded from the
claude-runtime; all actual HTTP calls are made by
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
   - buildOptions() with DATAGROK_PROMPT, tools, mcpServers, gate hooks
   - query({prompt, options}) → SDK event loop

5. Claude reads JS API source (Read tool on workspace/js-api/src/viewer.ts)
   → SDK fires assistant event (tool_use: Read ...)
   → forwardEvent() emits {type:'tool_activity', summary:'Read workspace/js-api/src/viewer.ts'}
   → browser onToolActivity → panel shows the badge

6. Claude calls datagrok_exec with the viewer code + inline verify assertion
   → runtime emits {type:'input_request', toolName:'datagrok_exec', input:{code, verify}}
   → browser executeSingleBlock() adds the scatter plot to the live view,
     runVerification() re-reads view.viewers, respondToInput({success:true, verified:{passed:true}})
   → the verify gate marks the action proven

7. Claude emits text delta chunks
   → stream_event → {type:'chunk', content:'Added a scatter plot…'}
   → browser onChunk → accumulated text, panel.updateStreaming()

8. SDK result event → {type:'final', content, metrics}
   → browser onFinal → panel.finalizeStreaming() renders the final markdown
```

---

## Key Files Reference

| File                                                                                      | Purpose                                                                        |
|-------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
| [`dockerfiles/claude-runtime/Dockerfile`](../dockerfiles/claude-runtime/Dockerfile)       | Container build: clones public repo to /workspace, compiles `src/` with tsc    |
| [`dockerfiles/claude-runtime/src/server.ts`](../dockerfiles/claude-runtime/src/server.ts) | Hono server: `/ws` dispatch, task claim/release, auth subprocess, provider config |
| [`dockerfiles/claude-runtime/src/session.ts`](../dockerfiles/claude-runtime/src/session.ts) | Session registry, turn queue, event forwarding, browser round-trips, abort   |
| [`dockerfiles/claude-runtime/src/query-options.ts`](../dockerfiles/claude-runtime/src/query-options.ts) | `buildOptions`, browser/view MCP tools, access guard, tool-activity summaries |
| [`dockerfiles/claude-runtime/src/prompts.ts`](../dockerfiles/claude-runtime/src/prompts.ts) | System prompts + `buildSystemPrompt`                                          |
| [`dockerfiles/claude-runtime/src/verify.ts`](../dockerfiles/claude-runtime/src/verify.ts) | Act → verify gate (`PostToolUse`/`Stop` hooks)                                 |
| [`dockerfiles/claude-runtime/src/grounding.ts`](../dockerfiles/claude-runtime/src/grounding.ts) | Grounding gate for platform how-to claims                                  |
| [`dockerfiles/claude-runtime/src/stream-filter.ts`](../dockerfiles/claude-runtime/src/stream-filter.ts) | Fence-aware chunk filtering, hidden revision streaming                 |
| [`dockerfiles/claude-runtime/src/tasks.ts`](../dockerfiles/claude-runtime/src/tasks.ts)   | Celery admission-task claim/release registry                                   |
| [`dockerfiles/claude-runtime/src/types.ts`](../dockerfiles/claude-runtime/src/types.ts)   | WebSocket message type definitions                                             |
| [`dockerfiles/mcp-server/src/index.ts`](../dockerfiles/mcp-server/src/index.ts)           | MCP tool registry                                                              |
| [`dockerfiles/mcp-server/src/api-client.ts`](../dockerfiles/mcp-server/src/api-client.ts) | HTTP client to Datagrok API                                                    |
| [`src/claude/runtime-client.ts`](../src/claude/runtime-client.ts)                         | Browser WebSocket client, RxJS subjects, `ClaudeRuntimeClient` singleton       |
| [`src/claude/exec-blocks.ts`](../src/claude/exec-blocks.ts)                               | `executeSingleBlock()`, `runVerification()`, `buildWorkspaceContext()`, entity cards |
| [`src/claude/queue-task.ts`](../src/claude/queue-task.ts)                                 | Queued-task admission (`aiChatTurnTask` holds the celery slot)                 |
| [`src/ai/view-tools.ts`](../src/ai/view-tools.ts)                                         | View-function meta-tools declared as `clientTools`                             |
| [`src/ai/ui.ts`](../src/ai/ui.ts)                                                         | `runPromptWithLifecycle()`, streaming orchestration, panel setup functions     |
| [`src/ai/panel.ts`](../src/ai/panel.ts)                                                   | `AIPanel` — the singleton chat panel (`StreamingPanel` interface)              |
