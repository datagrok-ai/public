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

## Architecture

### LLM Client Layer (`src/llm-utils/LLM-client.ts`)

`LLMClient` is the central singleton for all LLM communication. It wraps:

- A raw `OpenAI` client (used only for vector store search)
- `LanguageModelV3` instances from `@ai-sdk/*` providers, one per `ModelOption` (`Fast`, `Deep Research`, `Coding`)

The provider is chosen based on `PackageSettings.APIName` (enum `LLMApiNames`). Model names are configurable via package
settings and can be remapped by Azure deployment names.

All LLM requests go through the Datagrok AI proxy (`grok.ai.config.proxyUrl`), authenticated
with `grok.ai.config.proxyToken`.

### Credentials (`src/llm-utils/creds.ts`)

`LLMCredsManager` — singleton that holds the vector store ID from package settings. Initialized in the `init()`
function.

### AI Search Providers (`src/package.ts`)

Three `aiSearchProvider` functions are registered via decorators, each with a `useWhen` clause that
the `CombinedAISearchAssistant` uses to route user queries:

- **Help** — documentation Q&A via vector store file search
- **Execute** — plans and runs multi-step function chains (`ChatGptAssistant`)
- **Query** — finds and executes matching database queries

### Combined Search (`src/llm-utils/combined-search.ts`)

`CombinedAISearchAssistant` collects all `aiSearchProvider` functions, uses a fast LLM call to rank them by relevance to
the user query, then presents results in a lazy tab control.

### Prompt Engine (`src/prompt-engine/`)

- `PromptEngine` interface with two implementations: `ChatGPTPromptEngine` (cloud LLM) and `GeminiPromptEngine` (Chrome
  built-in `LanguageModel`)
- `ChatGptAssistant` — multi-step planner: selects relevant packages, discovers their functions, asks the LLM to produce
  a `Plan` (JSON schema-constrained), then executes steps sequentially resolving `$`-prefixed variable references
  between steps

### Table View AI (`src/llm-utils/tableview-tools.ts`)

Agentic tool-use loop for manipulating viewers in a `TableView`. Provides
tools: `add_viewer`, `adjust_viewer`, `describe_viewer`, `find_viewers_by_type`, `list_current_viewers`, `highlight_element`, `reply_to_user`, `search_documentation`.
Runs up to 15 iterations (extendable via user prompt).

### SQL Tools (`src/llm-utils/sql-tools.ts`)

Generates SQL queries from natural language using database schema metadata and tool calls.

### Claude Code Runtime (`src/claude-code/`)

The browser-facing layer for the Claude runtime.

- **`ClaudeRuntimeClient`** — singleton WebSocket client. Discovers `claude-runtime` and `mcp-server` containers via `grok.dapi.docker`, then exposes RxJS subjects for streaming events. API: `send()`, `abort()`, `respondToInput()`, and promise-based `query()` for one-shot calls.
- **`claude-panel.ts`** — processes two special fenced blocks in Claude responses: `datagrok-exec` (runs JS with `grok`/`ui`/`DG`/`view`/`t` globals) and `datagrok-entities` (renders interactive entity cards).

### AI Panels (`src/llm-utils/panel.ts`)

Three panel classes handle UI for different contexts:

- `TVAIPanel` — table view assistant (supports DeepGROK and Claude engines)
- `DBAIPanel` — database query editor assistant
- `ScriptingAIPanel` — script generation assistant

All panels share a common base with chat history, model selection, streaming display, and conversation persistence
via `ConversationStorage`.

### Docker Containers (`dockerfiles/`)

- `claude-runtime/` — Hono + WebSocket server wrapping `@anthropic-ai/claude-agent-sdk`. Runs Claude sessions with a Datagrok-specific system prompt and resumable session support. Streams structured events to the browser: `chunk`, `tool_activity`, `tool_result`, `final`, `error`, `aborted`, `input_request`.
- `mcp-server/` — MCP server (`@modelcontextprotocol/sdk`, HTTP transport) exposing Datagrok operations as tools: functions (list/get/call/create), files (list/download/upload), projects, spaces, and user info. Auth via per-request `x-user-api-key` / `x-datagrok-api-url` headers.

## Key Singletons

Most core classes use singleton pattern with `getInstance()`:

- `LLMClient.getInstance()`
- `LLMCredsManager.getInstance()`
- `CombinedAISearchAssistant.instance`
- `ChatGPTPromptEngine.getInstance()`
- `ClaudeRuntimeClient.getInstance()`

## Package Settings

Configured via Datagrok package properties (in `package.json`):

| Setting | Purpose |
|---------|---------|
| `APIName` | LLM provider: `OpenAI Chat Completions`, `OpenAI Responses`, `Anthropic Messages`, or `Amazon Bedrock` |
| `defaultFastModel` | Model for simple tasks (default: `gpt-4o-mini`) |
| `defaultDeepResearchModel` | Model for complex tasks (default: `gpt-5.2`) |
| `defaultCodingModel` | Model for code generation (default: `gpt-5.1-codex-max`) |
| `vectorStoreId` | OpenAI vector store ID for documentation search |
| `region` | AWS region (Bedrock only) |
| `apiVersion` | API version path segment |