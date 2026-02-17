# Claude Runtime Implementation Specification

## 1. Overview

Dockerized Claude Code runtime integrated into the Datagrok platform via the ChatGpt package.
Runs a Python Quart server inside a Docker container, communicates with the JS client over
WebSocket (proxied through Datagrok backend), and uses the Claude Agent SDK
(`claude-agent-sdk`) to execute agentic queries with built-in tools (file read/write/edit,
bash, grep, glob, web search, etc.) **plus custom Datagrok tools** (add/adjust viewers,
explore DB schemas, generate SQL, generate scripts) exposed via an MCP bridge that proxies
tool calls over WebSocket to the browser for execution. Streams responses and tool-use
events into the existing AIPanel UI.

## 2. File Structure

```
packages/ChatGpt/
├── dockerfiles/
│   └── claude-runtime/
│       ├── Dockerfile
│       ├── container.json
│       ├── requirements.txt
│       ├── app.py
│       └── datagrok_mcp.py          # MCP tool definitions for Datagrok browser tools
├── src/
│   ├── claude-code/
│   │   ├── claude-runtime-client.ts  # WebSocket client + browser-side tool executor
│   │   ├── claude-panel.ts           # AIPanel integration
│   │   └── tool-registry.ts          # Maps tool names → existing context class methods
│   └── llm-utils/
│       ├── tableview-tools.ts        # TableViewContext (reused as-is)
│       ├── sql-tools.ts              # SQLGenerationContext (reused as-is)
│       └── script-tools.ts           # ScriptGenerationContext (reused as-is)
```

Note: `claude_client.py` and `session_manager.py` from the original design are eliminated —
the Claude Agent SDK handles both the API client and session state management internally.

The existing `*-tools.ts` files contain context classes (`TableViewContext`,
`SQLGenerationContext`, `ScriptGenerationContext`) whose methods are **reused as-is** for
browser-side execution. The agent loops (`processAgentMode`, `generateAISqlQueryWithTools`,
etc.) in those files are not used by the Claude path — the Agent SDK runs its own tool loop.

## 3. Docker Container

### 3.1 Dockerfile

Based on the established pattern from `packages/CVMTests/dockerfiles/cvmtests-docker-test2/`.

- Base image: `python:3.11-slim` (multi-stage build; cannot use Alpine — the Agent SDK
  bundles a Node.js CLI binary that requires glibc)
- Stage 1 (`pip_builder`): create venv, install dependencies from `requirements.txt`
- Stage 1 also: install Node.js 20.x (required by `claude-agent-sdk` internally)
- Stage 2: copy venv and Node.js, create non-root `grok` user (UID 1001, GID 2001), copy app files
- Port: `5353`
- Healthcheck: `wget --spider http://127.0.0.1:5353/health` (interval 5s, timeout 3s, retries 3)
- Entrypoint: `python app.py`
- Env vars:
  - `ANTHROPIC_API_KEY` — Anthropic API key (required, passed at deploy time; this is the
    env var name the Agent SDK reads)
  - `PYTHONUNBUFFERED=1` — immediate log output

### 3.2 container.json

```json
{
  "cpu": 1,
  "memory": 3072,
  "on_demand": true,
  "gpu": 0
}
```

Memory bumped to 3 GiB — the Agent SDK spawns a Node.js subprocess per query which adds
overhead compared to a plain API client. Anthropic recommends ~1 GiB minimum per SDK process.

### 3.3 requirements.txt

```
quart==0.19.6
hypercorn==0.17.3
websockets==12.0
claude-agent-sdk>=1.0.0
```

Note: `anthropic` Python SDK is no longer needed — the Agent SDK replaces it entirely and
manages the API communication internally.

## 4. Python Application

### 4.1 app.py — Quart Server

Responsibilities:
- Initialize Quart app and Hypercorn config bound to `0.0.0.0:5353`
- Register error handlers (404, 500) returning JSON
- Expose endpoints:
  - `GET /health` — returns `{"status": "ok"}`, status 200
  - `WS /ws` — bidirectional WebSocket handler

WebSocket handler logic:
1. Receive JSON message from client
2. Parse `sessionId` and `message` fields
3. Create a `DatagrokMcpBridge` instance for this WebSocket connection (see §4.3)
4. Call `claude_agent_sdk.query()` with:
   - `prompt = message`
   - `options.resume = sessionId` (for follow-up messages in an existing conversation;
     omit on first message — the SDK returns a `session_id` in the init message, store the
     mapping `clientSessionId → sdkSessionId` in a dict)
   - `options.system_prompt` — Datagrok-specific system prompt
   - `options.allowed_tools` — SDK built-in tools + MCP tool names
     (e.g., `["Read", "Glob", "Grep", "WebSearch", "mcp__datagrok__add_viewer", ...]`)
   - `options.mcp_servers = {"datagrok": bridge.mcp_server}` — registers the Datagrok
     MCP server that proxies tool calls to the browser
   - `options.permission_mode = "bypassPermissions"` (server-side, no human in the loop)
   - `options.max_turns` — safety limit (e.g., 15)
   - `options.include_partial_messages = True` — enables streaming
5. Iterate over the async generator from `query()`:
   - On `SDKSystemMessage` (subtype `init`): store `sdkSessionId` mapping
   - On `SDKPartialAssistantMessage` (text deltas): send `{"type": "chunk", ...}`
   - On `SDKAssistantMessage` with tool-use content blocks: send `{"type": "tool_use", ...}`
   - On `SDKResultMessage` (subtype `success`): send `{"type": "final", ...}` with the
     last assistant message content, plus `usage` and `cost` metadata
   - On `SDKResultMessage` (error subtypes): send `{"type": "error", ...}`
6. On exception: send `{"type": "error", ...}`, log traceback

State:
- `dict[str, str]` mapping client `sessionId` → SDK `session_id` (for session resumption)
- Max 200 entries, evict oldest on overflow
- No need to track conversation history — the SDK manages it internally

Split into two files: `app.py` (~100 lines) and `datagrok_mcp.py` (~150 lines).

### 4.2 Agent SDK Integration

**Replaces:** `claude_client.py` (Anthropic streaming client) and `session_manager.py`
(in-memory session state). Both are fully handled by the Agent SDK.

#### Query invocation (inside WebSocket handler)

```python
from claude_agent_sdk import query, ClaudeAgentOptions
from datagrok_mcp import DatagrokMcpBridge

# Create MCP bridge for this WebSocket connection
bridge = DatagrokMcpBridge(ws)

options = ClaudeAgentOptions(
    system_prompt="You are a coding assistant for the Datagrok platform...",
    allowed_tools=[
        # SDK built-in tools
        "Read", "Glob", "Grep", "WebSearch", "WebFetch",
        # Datagrok MCP tools (prefixed with mcp__datagrok__)
        "mcp__datagrok__add_viewer",
        "mcp__datagrok__describe_viewer",
        "mcp__datagrok__list_current_viewers",
        "mcp__datagrok__adjust_viewer",
        "mcp__datagrok__try_sql",
        "mcp__datagrok__describe_tables",
        # ... etc.
    ],
    mcp_servers={"datagrok": bridge.mcp_server},
    permission_mode="bypassPermissions",
    max_turns=15,
    model="sonnet",
    include_partial_messages=True,
)

# First message — new session
async for message in query(prompt=user_text, options=options):
    ...  # forward to WebSocket

# Follow-up message — resume existing session
options.resume = sdk_session_id
async for message in query(prompt=user_text, options=options):
    ...  # forward to WebSocket
```

#### System prompt construction
- Static preamble: role description ("You are a coding assistant for the Datagrok platform...")
- Datagrok tool usage instructions: when to use `add_viewer` vs `try_sql`, etc.
  (adapted from existing system prompts in `tableview-tools.ts` and `sql-tools.ts`)
- Injected repo context: key file paths and summaries (from `context.repoFiles` in the
  incoming message, or loaded once at startup)
- Keep total system prompt under 8K tokens

#### Tool whitelist
Two categories of tools — SDK built-in and Datagrok MCP:

| Category | Tools |
|----------|-------|
| SDK built-in (read-only) | `Read`, `Glob`, `Grep`, `WebSearch`, `WebFetch` |
| SDK built-in (code assist) | Add `Edit`, `Write` |
| SDK built-in (full agent) | Add `Bash`, `Task`, `NotebookEdit` |
| Datagrok Table View | `mcp__datagrok__add_viewer`, `mcp__datagrok__describe_viewer`, `mcp__datagrok__adjust_viewer`, `mcp__datagrok__list_current_viewers`, `mcp__datagrok__find_viewers_by_type`, `mcp__datagrok__get_available_viewers`, `mcp__datagrok__highlight_element`, `mcp__datagrok__list_ui_elements` |
| Datagrok SQL | `mcp__datagrok__list_tables_in_schema`, `mcp__datagrok__describe_tables`, `mcp__datagrok__list_joins`, `mcp__datagrok__try_sql`, `mcp__datagrok__list_all_schemas`, `mcp__datagrok__list_catalogs` |
| Datagrok Scripts | `mcp__datagrok__find_similar_script_samples`, `mcp__datagrok__list_js_members`, `mcp__datagrok__search_documentation` |

The `allowed_tools` list is assembled dynamically based on the active view type and context
sent by the browser in the `context.metadata` field (e.g., `viewType: "TABLE_VIEW"` enables
table view tools; `connectionId: "..."` enables SQL tools).

#### Session management
The SDK manages conversation history internally. The server only maintains a lightweight
mapping `dict[str, str]` from client `sessionId` to SDK `session_id`:
- On first message: call `query()` without `resume`, capture `session_id` from the
  `SDKSystemMessage` (subtype `init`)
- On follow-up: pass `resume=sdk_session_id`
- On "New Chat": delete mapping entry, next message starts a fresh session
- Max 200 mapping entries, evict oldest on overflow
- No persistence — lost on container restart (acceptable for v1)

### 4.3 datagrok_mcp.py — MCP Bridge to Browser

This module defines an in-process MCP server whose tools are **proxied to the browser**
over WebSocket for execution. It bridges the gap between the Agent SDK (running in the
Docker container) and the Datagrok platform tools (running in the user's browser).

#### Architecture

```
┌─────────────────────────────────────┐     ┌──────────────────────────────────┐
│         Docker Container            │     │           Browser                │
│                                     │     │                                  │
│  Agent SDK query()                  │     │  ClaudeRuntimeClient             │
│    │                                │     │    │                             │
│    ├─ built-in: Read, Grep, Glob    │     │    │                             │
│    │  (executed locally)            │     │    │                             │
│    │                                │     │    │                             │
│    └─ MCP: datagrok-tools           │     │    ├─ tool_execute handler        │
│        │                            │     │    │   ├─ TableViewContext        │
│        ├─ add_viewer ──────────────────WS────→ │   ├─ SQLGenerationContext    │
│        ├─ describe_viewer ─────────────WS────→ │   └─ ScriptGenerationContext │
│        ├─ try_sql ─────────────────────WS────→ │                             │
│        ├─ list_current_viewers ────────WS────→ │                             │
│        └─ ... (all Datagrok tools) ────WS────→ │                             │
│                                     │     │    │                             │
│    streams chunks/final ───────────────WS────→ └─ updates AIPanel UI         │
└─────────────────────────────────────┘     └──────────────────────────────────┘
```

#### Class: `DatagrokMcpBridge`

Constructor takes the active WebSocket connection. Creates an MCP server with tool
definitions matching the existing tools from `tableview-tools.ts`, `sql-tools.ts`, and
`script-tools.ts`.

```python
class DatagrokMcpBridge:
    def __init__(self, ws):
        self.ws = ws
        self.pending: dict[str, asyncio.Future] = {}
        self.mcp_server = self._build_mcp_server()

    def _build_mcp_server(self):
        """Create MCP server with all Datagrok tool definitions."""
        tools = [
            # Table View tools (from tableview-tools.ts → TableViewContext)
            tool("add_viewer",
                 "Add a new viewer to the table view with specified properties. "
                 "Use describe_viewer first to understand available properties.",
                 {"viewerType": str, "viewerProperties": dict},
                 self._proxy),
            tool("describe_viewer",
                 "Get detailed info about a viewer type including all its properties.",
                 {"viewerType": str},
                 self._proxy),
            tool("list_current_viewers",
                 "List all currently added viewers with IDs, types, and properties.",
                 {},
                 self._proxy),
            tool("adjust_viewer",
                 "Adjust properties of an existing viewer by its ID.",
                 {"viewerId": str, "viewerProperties": dict},
                 self._proxy),
            tool("find_viewers_by_type",
                 "Find all viewers of a specific type and get their IDs and properties.",
                 {"viewerType": str},
                 self._proxy),
            tool("get_available_viewers",
                 "Get list of all available viewer types.",
                 {},
                 self._proxy),
            tool("highlight_element",
                 "Highlight a UI element to guide the user.",
                 {"elementId": str, "description": str},
                 self._proxy),
            tool("list_ui_elements",
                 "List all available UI elements with IDs and descriptions.",
                 {},
                 self._proxy),

            # SQL tools (from sql-tools.ts → SQLGenerationContext)
            tool("list_tables_in_schema",
                 "List all table names in a schema with descriptions.",
                 {"schemaName": str, "catalogName": str},
                 self._proxy),
            tool("describe_tables",
                 "Get detailed column info for tables. Use schema.table format.",
                 {"tables": list},
                 self._proxy),
            tool("list_joins",
                 "List foreign key relationships for specified tables.",
                 {"tables": list},
                 self._proxy),
            tool("try_sql",
                 "Execute SQL query to test it. Returns row count and columns (limited to 10 rows).",
                 {"sql": str, "description": str},
                 self._proxy),
            tool("list_all_schemas",
                 "List all schemas in a catalog.",
                 {"catalogName": str},
                 self._proxy),
            tool("list_catalogs",
                 "List all database catalogs available on the connection.",
                 {},
                 self._proxy),
            tool("find_similar_queries",
                 "Find similar previously executed queries based on a prompt.",
                 {"prompt": str},
                 self._proxy),

            # Script tools (from script-tools.ts → ScriptGenerationContext)
            tool("find_similar_script_samples",
                 "Find similar script code samples from Datagrok codebase.",
                 {"description": str},
                 self._proxy),
            tool("search_documentation",
                 "Search Datagrok documentation and API references.",
                 {"query": str, "maxResults": int},
                 self._proxy),
            tool("list_js_members",
                 "List members of a JS namespace to verify API existence.",
                 {"expression": str},
                 self._proxy),
            tool("search_js_api_sources",
                 "Search Datagrok JS API source code to verify a method exists.",
                 {"query": str, "maxResults": int},
                 self._proxy),
        ]
        return createSdkMcpServer(name="datagrok", tools=tools)

    async def _proxy(self, tool_name: str, args: dict) -> dict:
        """Send tool call to browser over WebSocket, await result."""
        call_id = str(uuid4())
        future = asyncio.get_event_loop().create_future()
        self.pending[call_id] = future

        await self.ws.send(json.dumps({
            "type": "tool_execute",
            "callId": call_id,
            "tool": tool_name,
            "args": args,
        }))

        try:
            result = await asyncio.wait_for(future, timeout=30.0)
        except asyncio.TimeoutError:
            result = f"Error: browser did not respond to {tool_name} within 30s"
        finally:
            self.pending.pop(call_id, None)

        return {"content": [{"type": "text", "text": result}]}

    def resolve(self, call_id: str, result: str):
        """Called when browser sends back a tool_result message."""
        future = self.pending.get(call_id)
        if future and not future.done():
            future.set_result(result)
```

#### Timeout and error handling
- Each browser tool call has a 30-second timeout (configurable)
- If the browser disconnects mid-tool-call, all pending futures are resolved with an error
- If the browser returns an error string, it is forwarded to the SDK as-is
- The SDK treats tool errors as informational (Claude can decide to retry or adjust)

## 5. WebSocket API Contract

### 5.1 Incoming Messages (JS → Container)

#### User message (initiates a query)
```json
{
  "type": "user_message",
  "sessionId": "string",
  "message": "string",
  "context": {
    "viewType": "TABLE_VIEW",
    "connectionId": "optional-db-connection-id",
    "repoFiles": ["string"],
    "metadata": {}
  }
}
```

- `sessionId` — unique per conversation, generated client-side (e.g. `DG.toDart()` UUID or `crypto.randomUUID()`)
- `message` — user's text input
- `context.viewType` — active Datagrok view type; determines which MCP tools are enabled
  (e.g., `TABLE_VIEW` enables viewer tools, presence of `connectionId` enables SQL tools)
- `context.connectionId` — if set, enables SQL tools scoped to this database connection
- `context.repoFiles` — optional list of file paths to inject into prompt
- `context.metadata` — optional key-value pairs for prompt augmentation

#### Tool result (response to a `tool_execute` request from the container)
```json
{
  "type": "tool_result",
  "callId": "uuid",
  "result": "string (tool output or error message)"
}
```

- `callId` — matches the `callId` from the `tool_execute` request
- `result` — string output from executing the tool in the browser (same format as the
  existing tool return values in `tableview-tools.ts`, `sql-tools.ts`, etc.)

### 5.2 Outgoing Messages (Container → JS)

Streaming chunk (sent per token batch, from `SDKPartialAssistantMessage` text deltas):
```json
{
  "type": "chunk",
  "sessionId": "string",
  "content": "partial text"
}
```

Tool use notification (sent when the SDK executes a built-in tool — informational only):
```json
{
  "type": "tool_use",
  "sessionId": "string",
  "tool": "Read",
  "input": {"file_path": "/src/package.ts"},
  "status": "running"
}
```

Tool execute request (sent when the SDK calls a Datagrok MCP tool — requires browser
to execute the tool and respond with `tool_result`):
```json
{
  "type": "tool_execute",
  "callId": "uuid",
  "tool": "add_viewer",
  "args": {"viewerType": "Scatter plot", "viewerProperties": {"xColumnName": "Age"}}
}
```

- `callId` — unique ID; the browser must include this in its `tool_result` response
- `tool` — Datagrok tool name (matches method names in `TableViewContext`, `SQLGenerationContext`, etc.)
- `args` — tool arguments (matches the existing tool input schemas)

The container **blocks** on the SDK's tool execution until the browser responds with a
matching `tool_result` (see §5.1). Timeout: 30 seconds.

Final message (sent once after the SDK query completes with `SDKResultMessage` subtype `success`):
```json
{
  "type": "final",
  "sessionId": "string",
  "content": "full accumulated response",
  "usage": {"input_tokens": 1234, "output_tokens": 567},
  "cost_usd": 0.042,
  "num_turns": 3
}
```

Error message (sent on `SDKResultMessage` with error subtype, or on exception):
```json
{
  "type": "error",
  "sessionId": "string",
  "message": "error description"
}
```

## 6. TypeScript Client

### 6.1 claude-runtime-client.ts

Location: `src/claude-code/claude-runtime-client.ts`

Class: `ClaudeRuntimeClient`

Imports: `grok`, `DG`, `rxjs`

#### Properties

```typescript
private ws: WebSocket | null;
private container: DG.DockerContainer | null;
private toolExecutor: BrowserToolExecutor | null;
public onChunk: rxjs.Subject<{sessionId: string, content: string}>;
public onToolUse: rxjs.Subject<{sessionId: string, tool: string, input: any, status: string}>;
public onFinal: rxjs.Subject<{sessionId: string, content: string, usage?: any, cost_usd?: number, num_turns?: number}>;
public onError: rxjs.Subject<{sessionId: string, message: string}>;
```

#### Methods

`async connect(): Promise<void>`
1. Find container: `grok.dapi.docker.dockerContainers.filter('name = "claude-runtime"').first()`
2. If container not found or not running, throw with clear error message
3. Open WebSocket: `grok.dapi.docker.dockerContainers.webSocketProxy(container.id, '/ws', 120000)` (120s timeout for on_demand cold start)
4. Create `BrowserToolExecutor` (see §6.3)
5. Attach `ws.onmessage` handler:
   - Parse JSON
   - Route by `type` field:
     - `chunk` → `onChunk.next(...)`
     - `tool_use` → `onToolUse.next(...)` (informational, for UI display)
     - `tool_execute` → `toolExecutor.execute(data)` then send `tool_result` back
     - `final` → `onFinal.next(...)`
     - `error` → `onError.next(...)`
6. Attach `ws.onclose` / `ws.onerror` — log, set `ws = null`

`send(sessionId: string, message: string, context?: {viewType?: string, connectionId?: string, repoFiles?: string[], metadata?: Record<string, any>}): void`
1. If `ws` is null or not OPEN, throw
2. `ws.send(JSON.stringify({type: 'user_message', sessionId, message, context: context ?? {}}))`

`disconnect(): void`
1. Close WebSocket if open
2. Complete all Subjects
3. Set `ws = null`

`get connected(): boolean`
- Return `ws !== null && ws.readyState === WebSocket.OPEN`

Singleton pattern: `static getInstance(): ClaudeRuntimeClient`

### 6.2 claude-panel.ts

(unchanged — see existing section above)

### 6.3 tool-registry.ts — Browser-Side Tool Executor

Location: `src/claude-code/tool-registry.ts`

Class: `BrowserToolExecutor`

This class receives `tool_execute` messages from the container and dispatches them to the
existing context classes that already implement all the tool logic. It is the browser-side
counterpart of `DatagrokMcpBridge` in the container.

#### Constructor

Takes the active `DG.TableView` (if any) and `connectionId` (if any). Creates context
instances on demand:

```typescript
class BrowserToolExecutor {
  private tableViewContext: TableViewContext | null;
  private sqlContext: SQLGenerationContext | null;
  private scriptContext: ScriptGenerationContext | null;

  constructor(
    private tableView: DG.TableView | null,
    private connectionId: string | null,
  ) {
    this.tableViewContext = tableView ? new TableViewContext(tableView) : null;
    // sqlContext and scriptContext created lazily when needed
  }
```

#### Method: `async execute(msg: {callId: string, tool: string, args: any}): Promise<string>`

Routes by tool name to the appropriate context class method. The mapping is direct —
tool names match method names:

```typescript
async execute(msg: {callId: string, tool: string, args: any}): Promise<string> {
  const {tool, args} = msg;
  switch (tool) {
    // Table View tools → TableViewContext
    case 'add_viewer':
      return this.tableViewContext!.addViewer(args.viewerType, args.viewerProperties);
    case 'describe_viewer':
      return this.tableViewContext!.describeViewer(args.viewerType);
    case 'list_current_viewers':
      return this.tableViewContext!.listCurrentViewers();
    case 'adjust_viewer':
      return this.tableViewContext!.adjustViewer(args.viewerId, args.viewerProperties);
    case 'find_viewers_by_type':
      return this.tableViewContext!.findViewersByType(args.viewerType);
    case 'get_available_viewers':
      return this.tableViewContext!.getAvailableViewers();
    case 'highlight_element':
      return (await this.tableViewContext!.highlightElement(args.elementId, args.description))
        ?? 'Element highlighted';
    case 'list_ui_elements':
      return this.tableViewContext!.listUIElements();

    // SQL tools → SQLGenerationContext
    case 'list_tables_in_schema':
      return (await this.getSqlContext().listTables(false, args.schemaName, args.catalogName)).description;
    case 'describe_tables':
      return this.getSqlContext().describeTables(args.tables);
    case 'list_joins':
      return this.getSqlContext().listJoins(args.tables);
    case 'try_sql':
      return this.getSqlContext().trySql(args.sql, args.description);
    case 'list_all_schemas':
      return (await this.getSqlContext().listAllSchemas(args.catalogName)).join(', ');
    case 'list_catalogs':
      return (await this.getSqlContext().listCatalogs()).join(', ');

    // Script tools → ScriptGenerationContext
    case 'find_similar_script_samples':
      return this.getScriptContext().findSimilarScriptSamples(args.description, args.language);
    case 'search_documentation':
      return this.getScriptContext().searchDocumentation(args.query, args.maxResults);
    case 'list_js_members':
      return this.getScriptContext().listJsMembers(args.expression);
    case 'search_js_api_sources':
      return this.getScriptContext().searchJsApiSources(args.query, args.maxResults);

    default:
      return `Error: unknown tool "${tool}"`;
  }
}
```

#### Integration in `ws.onmessage`

When a `tool_execute` message arrives:

```typescript
case 'tool_execute': {
  const {callId, tool, args} = data;
  let result: string;
  try {
    result = await this.toolExecutor!.execute({callId, tool, args});
  } catch (e: any) {
    result = `Error: ${e.message}`;
  }
  this.ws!.send(JSON.stringify({type: 'tool_result', callId, result}));
  // Also notify UI about the tool execution
  this.onToolUse.next({sessionId: data.sessionId, tool, input: args, status: 'completed'});
  break;
}
```

#### Key design points
- **No duplication**: All tool logic lives in the existing context classes. The executor
  is a thin routing layer (~50 lines).
- **Existing classes are reused as-is**: `TableViewContext`, `SQLGenerationContext`, and
  `ScriptGenerationContext` require no modifications. Their public methods already return
  `string` or `Promise<string>` — exactly what the MCP bridge expects.
- **Context classes need to be visible**: The classes in `tableview-tools.ts` and
  `sql-tools.ts` are currently not exported. Add `export` to `TableViewContext` and
  `SQLGenerationContext` (they are already structured as standalone classes).
- **Lazy initialization**: SQL and script contexts are created only when the first tool
  call of that type arrives, since they may require async setup (e.g., loading DB metadata).

### 6.4 claude-panel.ts

Location: `src/claude-code/claude-panel.ts`

Function: `setupClaudeAIPanelUI()`

This function hooks the Claude runtime into Datagrok views. It follows the exact same pattern
as `setupTableViewAIPanelUI()` and `setupScriptsAIPanelUI()` in `src/llm-utils/ui.ts`.

Implementation steps:

1. Check `grok.ai.config.configured` — return early if not configured
2. Subscribe to `grok.events.onViewAdded`:
   - For each relevant view (start with `TABLE_VIEW`), create an `AIPanel` instance
   - Alternatively, create a dedicated `ClaudeAIPanel` subclass that overrides `placeHolder`
     and hides the model selector (since model is fixed to Claude)
3. Subscribe to `panel.onRunRequest`:
   - Get or create `ClaudeRuntimeClient` singleton
   - If not connected, call `await client.connect()`
   - Generate a `sessionId` (use view ID or `crypto.randomUUID()`)
   - Determine context: `viewType` from current view, `connectionId` if applicable
   - Call `panel.startChatSession()` to get `{session, endSession}`
   - Add user message to UI: `session.addUserMessage(...)`
   - Send message: `client.send(sessionId, prompt, {viewType, connectionId})`
   - Subscribe to `client.onChunk` filtered by `sessionId`:
     - Accumulate content
     - Call `session.addAIMessage(...)` to update UI with streaming text
   - Subscribe to `client.onToolUse` filtered by `sessionId`:
     - Show tool activity in UI (e.g., "Adding Scatter plot viewer...")
   - Subscribe to `client.onFinal` filtered by `sessionId`:
     - Call `session.addAIMessage(...)` with final content and
       `messageOptions: {finalResult: content}` (enables feedback buttons)
     - Call `endSession()`
     - Unsubscribe from chunk/final/error/toolUse
   - Subscribe to `client.onError` filtered by `sessionId`:
     - Show error via `grok.shell.error(...)`
     - Call `endSession()`

#### Streaming UI Update Strategy

The existing `AIPanel.appendMessage()` creates a new DOM element per call. For streaming,
this would create one element per token — unusable. Two approaches:

**Option A (recommended):** Add a `updateLastAIMessage(content: string)` method to
`AIPanelFuncs`:
- On first chunk: call `addAIMessage()` to create the DOM element
- On subsequent chunks: find the last `.d4-ai-assistant-response-container` in the output
  area and replace its markdown content
- On final: replace with final content including feedback buttons

**Option B:** Accumulate all chunks in a buffer, only call `addAIMessage()` on `final`.
Show a typing indicator during streaming. Simpler but no incremental display.

Choose Option A for better UX. This requires a small addition to `AIPanel` or direct DOM
manipulation in `claude-panel.ts`.

## 7. Package Registration

### 7.1 package.ts Changes

Add to the ChatGpt package entry point (`src/package.ts`):

```typescript
//name: setupClaudeRuntime
//description: Initializes Claude Code runtime integration
//meta.autostartImmediate: true
export async function setupClaudeRuntime() {
  await setupClaudeAIPanelUI();
}
```

Or call `setupClaudeAIPanelUI()` from an existing `autostart` function if one exists.

### 7.2 Package Properties

Add to `package.json` properties array:

```json
{
  "name": "claudeRuntimeEnabled",
  "propertyType": "bool",
  "description": "Enable Claude Code runtime Docker container integration",
  "defaultValue": false,
  "nullable": false
}
```

Check this property in `setupClaudeAIPanelUI()` before activating.

## 8. Session Routing

- Each browser tab / Datagrok view gets its own client-side `sessionId`
- `sessionId` is generated client-side on first message and reused for the conversation
- The Docker container maintains a lightweight mapping `clientSessionId → sdkSessionId`
- On first message: `query()` is called without `resume`; the SDK returns a `session_id` in
  the `SDKSystemMessage` (subtype `init`), which is stored in the mapping
- On follow-up messages: `query()` is called with `resume=sdkSessionId` — the SDK restores
  full conversation context internally (no manual history management)
- Multiple concurrent sessions are supported (different users, different tabs)
- No cross-session state — sessions are fully isolated
- On "New Chat" button click: generate a new client `sessionId`, delete the old mapping entry;
  the SDK session remains in its internal storage until garbage-collected

## 9. Event Integration

No custom global event bus. Use `rxjs.Subject` internally within `ClaudeRuntimeClient`:
- `onChunk`, `onFinal`, `onError` — already defined as Subjects
- Other package components can subscribe to these if needed
- Follows the same pattern as `AIPanel._onRunRequest` and `_onClearChatRequest`

If other packages need to listen to Claude events in the future, expose an observable:
```typescript
// In package.ts
//name: claudeResponseStream
//output: observable stream
export function claudeResponseStream() {
  return ClaudeRuntimeClient.getInstance().onFinal.asObservable();
}
```

## 10. Error Handling

### Container side (Python)
- Check `SDKResultMessage.subtype` for error conditions:
  - `"error_max_turns"` — agent hit the turn limit; send error with "Max turns reached"
  - `"error_tool"` — a tool execution failed; include tool name and error details
  - `"error_authentication"` — API key invalid; send error with "Authentication failed"
  - Other error subtypes — forward the SDK's error message
- Catch JSON parse errors on incoming messages — send error with "Invalid message format"
- Catch all unexpected exceptions — log traceback, send generic error
- Use `max_turns` option (e.g., 10) as a safety limit to prevent runaway agent loops
- Optionally set `max_budget_usd` per query to cap API spend
- Never crash the WebSocket loop — always continue accepting new messages

### Client side (TypeScript)
- `connect()` failure — show `grok.shell.error('Claude runtime container not available')`
- WebSocket close during streaming — call `endSession()`, show error
- JSON parse failure on incoming message — log and ignore
- Timeout on connect (120s) — container may not be deployed; show actionable error

## 11. Testing Strategy

### Container tests (Python)
- Unit test session mapping: add/get/evict/clear of `clientSessionId → sdkSessionId` dict
- Unit test `DatagrokMcpBridge`: verify `tool_execute` is sent over WebSocket, verify
  `tool_result` resolves the pending future, verify 30s timeout
- Integration test: start Quart app, connect WebSocket, send message, verify
  chunk/tool_use/tool_execute/final responses
- Test session resumption: send two messages with the same `sessionId`, verify second
  uses `resume`
- Test `max_turns` enforcement: verify agent stops after limit and returns error

### Client tests (TypeScript)
- Register tests in `src/tests/claude-runtime-tests.ts`
- Test `ClaudeRuntimeClient.connect()` — requires running container
- Test send/receive message flow
- Test `BrowserToolExecutor` routing: mock context classes, verify each tool name
  dispatches to the correct method
- Test `tool_execute` → `tool_result` round-trip over WebSocket
- Test reconnection on WebSocket close
- Add to `package-test.ts` exports

### Manual testing
- Deploy container, open Datagrok table view, toggle AI panel, type a message
- Verify streaming chunks appear incrementally
- Ask Claude to "add a histogram for the Age column" — verify `tool_execute` is sent,
  viewer appears in the table view, and Claude's response references the result
- Ask Claude to "show me tables in the public schema" — verify SQL tools work via the bridge
- Verify "New Chat" resets the conversation
- Verify error display when container is stopped

## 12. Security Considerations

- `ANTHROPIC_API_KEY` is passed as an environment variable at deploy time, never baked into the image
- WebSocket traffic is proxied through `grok.dapi.docker.dockerContainers.webSocketProxy` — inherits Datagrok session authentication
- No direct network access to the container from outside the platform
- Session data is ephemeral (managed by SDK internally), no PII persistence
- Input sanitization: validate incoming JSON schema before processing
- **Tool whitelist (`allowed_tools`)** is the primary security boundary — controls both
  SDK built-in tools and Datagrok MCP tools. Start with read-only tools and expand as needed
- `permission_mode="bypassPermissions"` is used because the server has no interactive user —
  tool access is controlled entirely via `allowed_tools`
- Datagrok MCP tools execute **in the browser** under the user's session — they inherit the
  user's Datagrok permissions (e.g., `try_sql` respects DB connection access rights)
- The `tool_execute` → `tool_result` bridge only accepts messages from the container's
  WebSocket — no external injection possible (proxied through Datagrok backend auth)
- `max_turns` and `max_budget_usd` options prevent runaway agent loops and cost overruns
- Rate limiting: defer to Anthropic API rate limits; optionally add per-session rate limiting
  in the WebSocket handler

## 13. Deployment

1. Build and publish the ChatGpt package with `grok publish --release`
2. The platform auto-detects `dockerfiles/claude-runtime/` and builds the container
3. Set `ANTHROPIC_API_KEY` environment variable in the Datagrok Docker container configuration
4. Enable `claudeRuntimeEnabled` in ChatGpt package settings
5. Container starts on-demand on first WebSocket connection (due to `"on_demand": true`)
6. The container needs outbound HTTPS access to `api.anthropic.com` for the Agent SDK

## 14. Built-in Capabilities (from Agent SDK)

The following features that were listed as "future improvements" in the original design are
now available out of the box via the Agent SDK:

- **Tool use** — file read/write/edit, bash execution, glob, grep, web search/fetch
  (controlled via `allowed_tools`)
- **Multi-model** — pass `model="haiku"`, `"sonnet"`, or `"opus"` per query
- **Code execution** — built-in `Bash` tool (add to `allowed_tools` when ready)
- **Session management** — native `resume` / `fork_session` / `continue` options
- **MCP servers** — define custom tools via `mcp_servers` option (e.g., Datagrok-specific
  tools for querying datasets, running scripts, etc.)
- **Token tracking** — `SDKResultMessage` includes `usage` (input/output tokens),
  `total_cost_usd`, `duration_ms`, `num_turns`
- **Subagents** — define specialized agents via the `agents` option that the main agent can
  delegate to via the `Task` tool

## 15. Remaining Future Improvements

- **Persistent sessions**: Use SDK session files + persistent volume or Redis for cross-restart continuity
- **RAG context**: Vector-embed repository files, retrieve relevant snippets per query
- **Container scaling**: Multiple container replicas with session affinity for high concurrency
- **File upload**: Allow users to attach files from Datagrok file browser as context
- **Additional MCP tools**: Extend the Datagrok MCP bridge with more tools as needed
  (e.g., "run script", "create column from formula", "apply filter", "export to CSV")
- **Cost dashboard**: Aggregate `cost_usd` from final messages into a per-user usage dashboard
- **Tool use UI**: Show collapsible tool-use events in the AI panel (file reads, viewer
  additions, SQL test results, etc.)
- **Dynamic tool selection**: Automatically enable/disable MCP tools based on which Datagrok
  packages are installed (e.g., Chem tools only if Chem package is loaded)

## 16. Implementation Order

1. **Docker container (basic)** — `Dockerfile` (Python + Node.js), `container.json`,
   `requirements.txt`, `app.py` with SDK `query()` but **no MCP tools yet**
2. **Test container locally** — `docker build`, `docker run`, test with `wscat`;
   verify Agent SDK `query()` works inside the container (requires `ANTHROPIC_API_KEY`)
3. **TypeScript client (basic)** — `claude-runtime-client.ts` with chunk/final/error handling
4. **Panel integration** — `claude-panel.ts`, streaming UI updates
5. **End-to-end test (basic)** — deploy, verify plain text chat works through the full stack
6. **MCP bridge (container side)** — `datagrok_mcp.py` with `DatagrokMcpBridge` class,
   wire into `app.py` via `mcp_servers` option
7. **Browser tool executor** — `tool-registry.ts` with `BrowserToolExecutor` class,
   export `TableViewContext` and `SQLGenerationContext`, wire `tool_execute` handler into
   `claude-runtime-client.ts`
8. **End-to-end test (tools)** — test viewer manipulation, SQL exploration, script search
   via the MCP bridge
9. **Package registration** — update `package.ts`, add package property
10. **Polish** — error messages, reconnection logic, `allowed_tools` tuning, tool timeout
    handling, edge cases