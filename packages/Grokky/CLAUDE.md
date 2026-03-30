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

### Claude Runtime (`src/claude-code/`)

`ClaudeRuntimeClient` connects via WebSocket to a Docker container (`grokky-claude-runtime`) for streaming Claude
responses. Supports chunk/tool-activity/tool-result/final/error/abort/input-request events. A
companion `grokky-mcp-server` container exposes an MCP endpoint.

### AI Panels (`src/llm-utils/panel.ts`)

Three panel classes handle UI for different contexts:

- `TVAIPanel` — table view assistant (supports DeepGROK and Claude engines)
- `DBAIPanel` — database query editor assistant
- `ScriptingAIPanel` — script generation assistant

All panels share a common base with chat history, model selection, streaming display, and conversation persistence
via `ConversationStorage`.

### Docker Containers (`dockerfiles/`)

- `claude-runtime/` — Node.js WebSocket server hosting Claude AI sessions
- `mcp-server/` — MCP (Model Context Protocol) server for Datagrok tool access

## Key Singletons

Most core classes use singleton pattern with `getInstance()`:

- `LLMClient.getInstance()`
- `LLMCredsManager.getInstance()`
- `CombinedAISearchAssistant.instance`
- `ChatGPTPromptEngine.getInstance()`
- `ClaudeRuntimeClient.getInstance()`

## Package Settings

Configured via Datagrok package properties (in `package.json`):

- `APIName` — LLM
  provider: `"OpenAI Chat Completions"` | `"OpenAI Responses"` | `"Anthropic Messages"` | `"Amazon Bedrock"`
- `defaultFastModel`, `defaultDeepResearchModel`, `defaultCodingModel` — model names per tier
- `vectorStoreId` — OpenAI vector store for documentation search
- `region` — AWS region (Bedrock only)
- `apiVersion` — API version path segment