# Grokky

Grokky is the [AI for Datagrok](../../help/explore/ai/ai.md), enabling the following:
1. A generic AI JS API for core and other plugins to use LLMs
2. [Datagrok MCP server](dockerfiles/mcp-server/) so that Datagrok can be used by others
3. [Claude Code environment](dockerfiles/claude-runtime) for
   - writing scripts
   - doing agentic work on behalf of the user, augmented with custom skills
   and company-specific knowledge on processes, infrastruture, etc.
4. A general-purpose [AI console](...) to work with the platform
5. Supercharging your daily work:
   - Natural-language database queries
   - Automatic visualizations
   - Search across everything
   - Talk to your documents

Ideal for enterprises. You deploy Datagrok in your private cloud, and provide 
credentials for your enterprise Claude Code instance. Your data and tool usage
is safe, secure, governed, monitored and measured.

See also [Grokky GitHub issues](https://github.com/datagrok-ai/public/issues/3710). 

## Source Structure

```
src/
├── package.ts              # Entry point — registers functions, init, search providers
├── utils.ts                # Shared utilities (viewer/dataframe descriptions, events)
├── ai/                     # AI panels, search, and UI wiring
│   ├── panel.ts            # TVAIPanel, DBAIPanel, ScriptingAIPanel, StreamingPanel
│   ├── ui.ts               # Setup functions that wire panels into the platform UI
│   ├── storage.ts          # ConversationStorage — IndexedDB persistence for chat history
│   ├── usage-limiter.ts    # Per-user daily request limits (group-configurable)
│   └── search/
│       ├── combined-search.ts   # Routes user queries to ranked aiSearchProvider functions
│       └── query-matching.ts    # LLM-based matching of natural-language to DB queries
├── claude/                 # Browser-facing Claude runtime integration
│   ├── runtime-client.ts   # WebSocket client to claude-runtime container
│   └── exec-blocks.ts      # Executes datagrok-exec / datagrok-entities fenced blocks
├── db/                     # Database tooling
│   ├── sql-tools.ts        # SQLGenerationContext — tool-call-based SQL generation
│   ├── db-index-tools.ts   # DBSchemaInfo caching, metadata generation
│   ├── query-meta-utils.ts # BuiltinDBInfoMeta — extends DG.DbInfo with built-in relations
│   └── indexes/            # Hardcoded schema indexes for known databases
│       ├── biologics-index.ts
│       └── chembl-index.ts
└── depr/                   # Deprecated code kept for reference (not imported by active code)
```

## Architecture

### AI Search Providers (`src/package.ts`)

Three `aiSearchProvider` functions are registered via decorators, each with a `useWhen` clause that
`CombinedAISearchAssistant` uses to route user queries:

- **Help** — documentation Q&A
- **Execute** — plans and executes JS code using structured output
- **Query** — matches natural-language to registered database queries

All three route through `ClaudeRuntimeClient` (Help uses plain `query()`, Execute and Query use `query()` with
a structured output schema).

### Combined Search (`src/ai/search/combined-search.ts`)

`CombinedAISearchAssistant` collects all `aiSearchProvider` functions, uses `ClaudeRuntimeClient.query()` to rank them
by relevance to the user query, then presents results in a lazy tab control.

### AI Panels (`src/ai/panel.ts`)

Three panel classes handle UI for different contexts:

- `TVAIPanel` — table view assistant
- `DBAIPanel` — database query editor assistant
- `ScriptingAIPanel` — script generation assistant

All panels stream through `ClaudeRuntimeClient` and share a common `StreamingPanel` base with chat history,
streaming display, and conversation persistence via `ConversationStorage`.

### AI UI Wiring (`src/ai/ui.ts`)

Setup functions called from `init()` that attach AI panels to platform UI elements:

- `setupSearchUI()` — wires `CombinedAISearchAssistant` into the global search bar
- `setupTableViewAIPanelUI()` — adds `TVAIPanel` to table views
- `setupScriptsAIPanelUI()` — adds `ScriptingAIPanel` to script views
- `setupAIQueryEditorUI()` — adds `DBAIPanel` to the query editor

Also exports `askWiki()` and `smartExecution()`, which back the Help and Execute search providers respectively.

### Usage Limiter (`src/ai/usage-limiter.ts`)

`UsageLimiter` — singleton that enforces per-user daily AI request limits. Limits are group-configurable
(e.g., "AI Folks" group gets 1000/day, default is 20/day). Counts are persisted in `grok.userSettings`.

### Conversation Storage (`src/ai/storage.ts`)

`ConversationStorage` — IndexedDB-backed persistence for chat conversations. Stores up to 20 conversations
per context (connection ID, view, etc.), each with full message history and UI state.

### Claude Runtime Client (`src/claude/runtime-client.ts`)

Singleton WebSocket client. Discovers `claude-runtime` and `mcp-server` containers via `grok.dapi.docker`,
then exposes RxJS subjects for streaming events. API: `send()`, `abort()`, `respondToInput()`, and
promise-based `query()` for one-shot structured-output calls.

### Executable Blocks (`src/claude/exec-blocks.ts`)

Processes two types of fenced blocks in Claude responses:
- `` ```datagrok-exec `` — runs JS with `grok`/`ui`/`DG`/`view`/`t` globals
- `` ```datagrok-entities `` — renders interactive entity cards (files, scripts, queries, connections, projects, spaces)

Also provides `buildViewContext()` — serializes the current view (table columns, script code) for Claude prompts.

### SQL Tools (`src/db/sql-tools.ts`)

`SQLGenerationContext` — manages tool execution for natural-language-to-SQL generation. Supports multiple
catalogs. Tool calls (`db_list_schemas`, `db_get_table_info`, etc.) are defined in the MCP server and forwarded from the Claude runtime to the browser for execution.

### DB Index Tools (`src/db/db-index-tools.ts`)

`DBSchemaInfo` — caches database schema metadata (table/column info, foreign key references) per connection+schema.
Also provides `genDBConnectionMeta()` for generating and persisting LLM-friendly schema descriptions.

### Docker Containers (`dockerfiles/`)

- `claude-runtime/` — Hono + WebSocket server wrapping `@anthropic-ai/claude-agent-sdk`. Runs Claude sessions with a Datagrok-specific system prompt and resumable session support. Streams structured events to the browser: `chunk`, `tool_activity`, `tool_result`, `final`, `error`, `aborted`, `input_request`.
- `mcp-server/` — MCP server (`@modelcontextprotocol/sdk`, HTTP transport) exposing Datagrok operations as tools: functions (list/get/call/create), files (list/download/upload), projects, spaces, and user info. Auth via per-request `x-user-api-key` / `x-datagrok-api-url` headers.

### Deprecated Code (`src/depr/`)

Contains the previous LLM abstraction layer (`LLMClient`, `ChatGptAssistant`, `PromptEngine`, etc.) and
earlier implementations of tableview tools, script tools, and AI panels. Kept for reference but not imported
by active code. All active AI features now route through `ClaudeRuntimeClient`.

## Key Singletons

- `ClaudeRuntimeClient.getInstance()` — WebSocket connection to Claude runtime
- `CombinedAISearchAssistant.instance` — AI search routing
- `UsageLimiter.getInstance()` — daily request limit enforcement

## Configuration

All LLM requests go through `ClaudeRuntimeClient`, which connects via WebSocket to the `claude-runtime`
Docker container. The container handles Claude API authentication internally.