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

- Base image: `node:20-slim`
- Clones the public Datagrok repo into `/workspace` at build time — this gives Claude read access to
  the full JS API source, packages, samples, and documentation without needing to fetch them at runtime.
- Runs `dist/server.js` (compiled from [`src/server.ts`](../dockerfiles/claude-runtime/src/server.ts))
  on port 5355.
- Requires `ANTHROPIC_API_KEY` in the environment.
- Optional: `MILVUS_TOKEN` / `OPENAI_API_KEY` to enable the `claude-context-mcp` vector-search sidecar.
- Health check: `GET /health` → `{status: "ok"}`.

**`mcp-server`** ([`Dockerfile`](../dockerfiles/mcp-server/Dockerfile)):

- Runs `dist/index.js` (compiled from [`src/index.ts`](../dockerfiles/mcp-server/src/index.ts))
  on port 3003.
- Exposes Datagrok operations as MCP tools (HTTP/streamable transport).
- Auth is passed per-request via `x-user-api-key` and `x-datagrok-api-url` headers forwarded from
  the claude-runtime on every call.

### Container Discovery (browser side)

[`ClaudeRuntimeClient.connect()`](../src/claude/runtime-client.ts#L44) discovers both containers at
runtime via `grok.dapi.docker`:

```typescript
// runtime-client.ts:49-58
const [runtimeContainers, mcpContainers] = await Promise.all([
  grok.dapi.docker.dockerContainers.filter('name = "grokky-claude-runtime"').list(),
  grok.dapi.docker.dockerContainers.filter('name = "grokky-mcp-server"').list(),
]);
this.ws = await grok.dapi.docker.dockerContainers.webSocketProxy(containerId, '/ws');
this.mcpServerUrl = `${grok.dapi.root}/docker/containers/proxy/${mcpContainerId}/mcp`;
```

The WebSocket and MCP HTTP traffic both flow through Datagrok's reverse proxy, so the browser never
talks to the containers directly.

---

## Message Protocol

All communication between the browser and `claude-runtime` is JSON over a single WebSocket connection.
Message types are defined in [`dockerfiles/claude-runtime/src/types.ts`](../dockerfiles/claude-runtime/src/types.ts).

### Incoming messages (browser → runtime)

| `type` | Purpose | Key fields |
|--------|---------|------------|
| `user_message` | Send a prompt, start or resume a session | `sessionId`, `message`, `apiKey`, `mcpServerUrl`, `outputSchema?` |
| `abort` | Cancel an in-progress query | `sessionId` |
| `input_response` | Reply to an `input_request` (tool call result or user answer) | `sessionId`, `value` |

### Outgoing messages (runtime → browser)

| `type` | Purpose | Key fields |
|--------|---------|------------|
| `chunk` | Incremental text delta from Claude | `sessionId`, `content` |
| `tool_activity` | Human-readable summary of a tool call in progress | `sessionId`, `summary` |
| `tool_result` | Result text returned by a tool | `sessionId`, `content` |
| `final` | Completed response (full text + optional structured output) | `sessionId`, `content`, `structured_output?` |
| `error` | Session error | `sessionId`, `message` |
| `aborted` | Confirmed abort | `sessionId` |
| `input_request` | Runtime needs the browser to execute a tool or answer a question | `sessionId`, `toolName`, `input` |

---

## Server-Side: claude-runtime

[`server.ts`](../dockerfiles/claude-runtime/src/server.ts) is a [Hono](https://hono.dev) server with a
single WebSocket endpoint at `/ws`.

### Session lifecycle

[`handleMessage()`](../dockerfiles/claude-runtime/src/server.ts#L286) handles each `user_message`:

1. Looks up an existing Claude session ID (stored in the in-memory `sessions` LRU map, max 200 entries)
   for continuation. If found, passes it as `resume` to the SDK so conversation history is preserved.
2. Calls `query()` from `@anthropic-ai/claude-agent-sdk` with:
   - The [system prompt](#system-prompt) (`DATAGROK_PROMPT`).
   - `allowedTools`: `Read`, `Glob`, `Grep`, `Edit`, `Write`, `Bash`, `WebSearch`, `WebFetch`,
     `AskUserQuestion`.
   - `mcpServers`: the `datagrok` MCP server (and optionally `claude-context` for vector search).
   - `permissionMode: 'acceptEdits'` — no confirmation prompts for file edits.
   - `model: 'opus'`.
3. Iterates the async event stream from the SDK, forwarding each event to the browser via
   [`forwardEvent()`](../dockerfiles/claude-runtime/src/server.ts#L231).

### Tool interception: `canUseTool`

[`canUseTool`](../dockerfiles/claude-runtime/src/server.ts#L307) is the hook that intercepts two
categories of tool calls before the SDK executes them:

- **`AskUserQuestion`** — sends an `input_request` to the browser, suspends the Claude query, and waits
  for an `input_response` with the user's answers.
- **`mcp__datagrok__db_*`** tools (`db_list_catalogs`, `db_list_schemas`, `db_list_tables`,
  `db_describe_tables`, `db_list_joins`, `db_try_sql`) — forwarded to the browser for client-side
  execution against the live database connection (see [SQL generation](#sql-generation)).

Both cases use the same async handshake: `input_request` out, `input_response` in, resolved by
[`handleInputResponse()`](../dockerfiles/claude-runtime/src/server.ts#L277).

### SDK event → WebSocket message mapping

[`forwardEvent()`](../dockerfiles/claude-runtime/src/server.ts#L231) maps SDK events:

| SDK event type               | → WS message      |
|------------------------------|-------------------|
| `system` (init)              | stores session ID |
| `stream_event` (text_delta)  | `chunk`           |
| `assistant` (tool_use block) | `tool_activity`   |
| `user` (tool result)         | `tool_result`     |
| `tool_progress`              | `tool_activity`   |
| `tool_use_summary`           | `tool_activity`   |
| `result` (success)           | `final`           |
| `result` (error)             | `error`           |

[`toolSummary()`](../dockerfiles/claude-runtime/src/server.ts#L200) generates human-readable labels for
tool activities using the formatters in `toolFormatters` and `mcpFormatters`.

### System prompt

[`DATAGROK_PROMPT`](../dockerfiles/claude-runtime/src/server.ts#L12) (lines 12–85) instructs Claude to:

- Always search the codebase at `$CLAUDE_WORKSPACE` before writing any API calls.
- Execute code exclusively via `` ```datagrok-exec `` fenced blocks (never plain `` ```javascript ``).
- Use `AskUserQuestion` to clarify ambiguous requests before acting.
- Reference Datagrok entities (files, scripts, queries, connections, projects, spaces) via
  `` ```datagrok-entities `` blocks so they render as interactive cards.
- Use available globals: `grok`, `ui`, `DG` (always), `view` (always), `t` (TableView only).

---

## Browser Side: ClaudeRuntimeClient

[`ClaudeRuntimeClient`](../src/claude/runtime-client.ts) is a singleton WebSocket client.

**API:**

| Method | Purpose |
|--------|---------|
| `send(sessionId, message, options?)` | Send a prompt; optionally request structured output via `outputSchema` |
| `abort(sessionId)` | Cancel an in-progress query |
| `respondToInput(sessionId, value)` | Reply to an `input_request` |
| `query(message, options?)` | One-shot promise-based call; wraps send/subscribe/cleanup |

**RxJS subjects** (one per outgoing message type):

```typescript
// runtime-client.ts:18-25
onChunk: Subject<ChunkEvent>
onToolActivity: Subject<ToolActivityEvent>
onToolResult: Subject<ToolResultEvent>
onFinal: Subject<FinalEvent>
onError: Subject<ErrorEvent>
onAborted: Subject<AbortedEvent>
onInputRequest: Subject<InputRequestEvent>
```

---

## Streaming to the Panel UI

[`runClaudeStreaming()`](../src/ai/ui.ts#L189) in [`src/ai/ui.ts`](../src/ai/ui.ts) is the central
orchestrator that connects `ClaudeRuntimeClient` events to the chat panel UI:

```
client.onChunk       → panel.updateStreaming(accumulated)          (live markdown update)
client.onToolActivity → panel.updateStreaming(accumulated + status) (tool badge)
client.onToolResult  → panel.updateStreaming(accumulated + result)  (collapsible result)
client.onFinal       → panel.finalizeStreaming(content, view)       (execute blocks, render cards)
client.onError       → grok.shell.error(msg) + cleanup
client.onAborted     → panel.clearStreaming() + cleanup
client.onInputRequest → clientToolHandler(toolName, input)          (DB tools)
                      → panel.showInputRequest(input)               (AskUserQuestion)
```

The `clientToolHandler` parameter is supplied by callers that need client-side tool execution (e.g.
[`setupAIQueryEditorUI()`](../src/ai/ui.ts#L141) passes `sqlContext.handleToolCall`). For all other
panels the default path handles `AskUserQuestion` by showing an inline form.

### Panel entry points

| Function | Panel class | Context |
|----------|-------------|---------|
| [`setupTableViewAIPanelUI()`](../src/ai/ui.ts#L299) | `TVAIPanel` | TableView — passes table name, columns, row count |
| [`setupAIQueryEditorUI()`](../src/ai/ui.ts#L141) | `DBAIPanel` | Query editor — passes connection ID, catalog names |
| [`setupScriptsAIPanelUI()`](../src/ai/ui.ts#L329) | `ScriptingAIPanel` | Script editor — passes current script code |
| [`setupSearchUI()`](../src/ai/ui.ts#L74) | (no panel) | Global search bar — routes to `CombinedAISearchAssistant` |

All entry points are called from `init()` in [`src/package.ts`](../src/package.ts).

---

## Converting Claude Responses to UI

Claude responses may contain two special fenced block types, processed by
[`src/claude/exec-blocks.ts`](../src/claude/exec-blocks.ts) after the `final` event.

### `datagrok-exec` blocks

```typescript
// exec-blocks.ts:5-23
export async function executeDatagrokBlocks(content: string, view: DG.ViewBase): Promise<HTMLElement[]>
```

Finds all `` ```datagrok-exec\n...\n``` `` blocks with a regex, then for each block:

1. Extracts the current `DG.DataFrame` from the view (if it is a `TableView`).
2. Executes the code as an async IIFE via `new Function(...)`, passing `grok`, `ui`, `DG`, `view`, `t`
   as arguments.
3. If the function returns an `HTMLElement`, appends it to the chat as a result widget.
4. On error: calls `grok.shell.error()` and continues to the next block.

Regular `` ```javascript `` blocks are rendered as code snippets only — they do not execute.

### `datagrok-entities` blocks

```typescript
// exec-blocks.ts:84-102
export function renderEntityBlocks(container: HTMLElement): void
```

Called after the markdown is rendered to the DOM. Finds all `<code class="language-datagrok-entities">`
elements, parses their text content as a JSON array of `DgEntityRef` objects, and replaces each `<pre>`
with a row of entity cards:

1. Renders a `ui.divText(ref.name)` placeholder immediately.
2. Fetches the actual entity from `grok.dapi.*` asynchronously.
3. On success, replaces the placeholder with `ui.renderCard(entity)` — an interactive card the user
   can click to open the entity.

Supported entity types: `file`, `script`, `query`, `connection`, `project`, `space`.

### `buildViewContext()`

[`buildViewContext()`](../src/claude/exec-blocks.ts#L25) serializes the current view state to prepend to
every prompt, giving Claude accurate context:

- **TableView**: `Table "name" (N rows): col1(type), col2(type), ...`
- **ScriptView**: full script source in a fenced block.

---

## SQL Generation

For natural-language database queries ([`DBAIPanel`](../src/ai/panel.ts)), DB schema tools are defined in
the MCP server but executed client-side in the browser, where the active database connection lives.

The runtime's `canUseTool` intercepts any `mcp__datagrok__db_*` call and sends an `input_request` to
the browser instead of letting the MCP server handle it. The browser's
[`runClaudeStreaming()`](../src/ai/ui.ts#L259) strips the `mcp__datagrok__` prefix and dispatches to
[`SQLGenerationContext.handleToolCall()`](../src/db/sql-tools.ts):

| Tool | Browser action |
|------|---------------|
| `db_list_catalogs` | Returns catalog names from connection parameters |
| `db_list_schemas` | Queries schema list from the active connection |
| `db_list_tables` | Queries table list for a schema |
| `db_describe_tables` | Returns detailed column/type metadata |
| `db_list_joins` | Returns foreign-key relationships |
| `db_try_sql` | Runs the SQL (with a `LIMIT 10`), returning column names and sample rows |

`db_try_sql` includes a safety check against destructive patterns (`DROP`, `DELETE`, `UPDATE`, `INSERT`,
`ALTER`, etc.) and shows a user confirmation dialog if any are detected.

[`SQLGenerationContext.create()`](../src/db/sql-tools.ts#L30) is called lazily on the first prompt in
a session — it fetches the connection and pre-loads all catalog metadata via `BuiltinDBInfoMeta`.

---

## MCP Server Tools

[`dockerfiles/mcp-server/src/index.ts`](../dockerfiles/mcp-server/src/index.ts) registers tools using
`@modelcontextprotocol/sdk`. Tool calls are authenticated per-request via headers forwarded from the
claude-runtime; all actual HTTP calls are made by
[`src/api-client.ts`](../dockerfiles/mcp-server/src/api-client.ts).

Tool categories:

| Category | Tools |
|----------|-------|
| Functions | `list_functions`, `get_function`, `call_function`, `list_scripts`, `list_queries`, `create_script`, `create_query` |
| Files | `list_files`, `download_file`, `upload_file` |
| Projects | `list_projects`, `get_project`, `create_project`, `delete_project`, … |
| Spaces | `list_spaces`, `get_space`, `create_space`, `create_subspace`, `list_space_children`, `read_space_file`, `write_space_file`, … |
| DB (client-side) | `db_list_catalogs`, `db_list_schemas`, `db_list_tables`, `db_describe_tables`, `db_list_joins`, `db_try_sql` |
| User | `get_current_user` |

---

## End-to-End Flow Example

Below is the full path for "Add a scatter plot of Height vs Weight":

```
1. User types in TVAIPanel → panel.onRunRequest fires
2. runPromptWithLifecycle() → UsageLimiter.tryCheckAndIncrement()
3. runClaudeStreaming() → client.ensureConnected()
4. ClaudeRuntimeClient.send(sessionId, prompt + table context)
   WS → {type:'user_message', sessionId, message, apiKey, mcpServerUrl}

5. claude-runtime handleMessage():
   - buildOptions() with DATAGROK_PROMPT, allowedTools, mcpServers
   - query({prompt, options}) → SDK event loop

6. Claude reads JS API source (Read tool on /workspace/js-api/src/viewer.ts)
   → SDK fires assistant event (tool_use: Read ...)
   → forwardEvent() emits {type:'tool_activity', summary:'Read js-api/src/viewer.ts'}
   → browser onToolActivity → panel shows "Read js-api/src/viewer.ts"

7. Claude emits text delta chunks
   → stream_event → {type:'chunk', content:'Here is...'}
   → browser onChunk → accumulated text, panel.updateStreaming()

8. Claude emits final response with datagrok-exec block
   → result event → {type:'final', content:'...\n```datagrok-exec\n...```'}
   → browser onFinal → panel.finalizeStreaming()

9. finalizeStreaming() calls:
   a. executeDatagrokBlocks() → runs JS, DG.Viewer.scatterPlot(t, {...}).root returned
      → HTMLElement appended to chat
   b. renderEntityBlocks() → no datagrok-entities in this response, no-op

10. Scatter plot appears inline in the chat panel.
```

---

## Key Files Reference

| File                                                                                      | Purpose                                                                        |
|-------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
| [`dockerfiles/claude-runtime/Dockerfile`](../dockerfiles/claude-runtime/Dockerfile)       | Container build: Node 20, clones public repo to /workspace                     |
| [`dockerfiles/claude-runtime/src/server.ts`](../dockerfiles/claude-runtime/src/server.ts) | Hono WebSocket server, system prompt, session management, event forwarding     |
| [`dockerfiles/claude-runtime/src/types.ts`](../dockerfiles/claude-runtime/src/types.ts)   | WebSocket message type definitions                                             |
| [`dockerfiles/mcp-server/src/index.ts`](../dockerfiles/mcp-server/src/index.ts)           | MCP tool registry                                                              |
| [`dockerfiles/mcp-server/src/api-client.ts`](../dockerfiles/mcp-server/src/api-client.ts) | HTTP client to Datagrok API                                                    |
| [`src/claude/runtime-client.ts`](../src/claude/runtime-client.ts)                         | Browser WebSocket client, RxJS subjects, `ClaudeRuntimeClient` singleton       |
| [`src/claude/exec-blocks.ts`](../src/claude/exec-blocks.ts)                               | `executeDatagrokBlocks()`, `renderEntityBlocks()`, `buildViewContext()`        |
| [`src/ai/ui.ts`](../src/ai/ui.ts)                                                         | `runClaudeStreaming()`, panel setup functions, `askWiki()`, `smartExecution()` |
| [`src/ai/panel.ts`](../src/ai/panel.ts)                                                   | `TVAIPanel`, `DBAIPanel`, `ScriptingAIPanel`, `StreamingPanel`                 |
| [`src/db/sql-tools.ts`](../src/db/sql-tools.ts)                                           | `SQLGenerationContext` — client-side DB tool execution                         |
