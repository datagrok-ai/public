"""Every Hermes API the driver uses, typed as the driver uses it. Hermes is imported on first use, so the unit
tests run without it."""
from collections.abc import Callable
from typing import Any, Literal, NotRequired, Protocol, TypedDict

from .protocol import JsonObject

ToolHandler = Callable[..., str]
Hook = Callable[..., Any]


class ToolSpec(TypedDict):
    name: str
    description: str
    parameters: JsonObject


class Connection(TypedDict):
    provider: str
    api_mode: str
    base_url: str
    # On the Anthropic Messages API Hermes sends a callable's key as Bearer, a string as x-api-key.
    api_key: str | Callable[[], str]


class ConversationResult(TypedDict, total=False):
    messages: list[JsonObject]
    final_response: str | None
    failed: bool
    error: str | None
    turn_exit_reason: str | None


class ApiError(TypedDict):
    type: str
    message: str


class Activity(TypedDict, total=False):
    last_activity_at: float | None
    last_activity_description: str


class ToolVerdict(TypedDict):
    action: Literal['block', 'modify']
    message: NotRequired[str]
    args: NotRequired[JsonObject]


class RequestOverride(TypedDict):
    request: JsonObject
    source: str


class Agent(Protocol):
    session_id: str
    max_iterations: int
    reasoning_config: JsonObject
    enabled_toolsets: list[str]
    tools: list[JsonObject]
    valid_tool_names: set[str]

    def run_conversation(self, user_message: str | list[JsonObject], conversation_history: list[JsonObject],
                         task_id: str) -> ConversationResult: ...

    def interrupt(self, message: str | None = None, *, hard_cancel: bool = False) -> bool: ...

    def clear_interrupt(self) -> bool: ...

    def get_activity_summary(self) -> Activity: ...


class SessionStore(Protocol):
    def get_messages_as_conversation(self, session_id: str, *, repair_alternation: bool = False) \
        -> list[JsonObject]: ...


class ToolRegistry(Protocol):
    def register(self, name: str, toolset: str, schema: ToolSpec, handler: ToolHandler) -> None: ...

    def deregister(self, name: str) -> None: ...

    def get_entry(self, name: str) -> object | None: ...


class PluginContext(Protocol):
    def register_hook(self, hook_name: str, callback: Hook) -> object: ...

    def register_middleware(self, kind: str, callback: Hook) -> object: ...


def new_agent(model: str, connection: Connection, session_id: str, store: SessionStore, toolsets: list[str],
              on_delta: Callable[[str | None], None]) -> Agent:
    from run_agent import AIAgent
    agent: Agent = AIAgent(model=model, **connection, quiet_mode=True, session_id=session_id, session_db=store,
                           enabled_toolsets=toolsets, skip_context_files=True, skip_memory=True,
                           skip_background_review=True, stream_delta_callback=on_delta)
    return agent


def open_store() -> SessionStore:
    from hermes_constants import get_hermes_home
    from hermes_state_registry import acquire
    store: SessionStore = acquire(get_hermes_home() / 'state.db')
    return store


def tool_registry() -> ToolRegistry:
    from tools.registry import registry
    tools: ToolRegistry = registry
    return tools


def tool_definitions(toolsets: list[str]) -> list[JsonObject]:
    from model_tools import get_tool_definitions
    definitions: list[JsonObject] = get_tool_definitions(enabled_toolsets=toolsets, quiet_mode=True)
    return definitions


def spillover_dir() -> str:
    from tools.tool_result_storage import get_spillover_dir
    return str(get_spillover_dir())


def set_terminal_cwd(task_id: str, cwd: str) -> None:
    from tools.terminal_tool import register_task_env_overrides
    register_task_env_overrides(task_id, {'cwd': cwd})


def connect_mcp(name: str, config: JsonObject) -> list[str]:
    """Every MCP tool Hermes now has: the server's are missing when it could not connect."""
    from tools.mcp_tool_discovery import register_mcp_servers
    names: list[str] = register_mcp_servers({name: config})
    return names


def disconnect_mcp(name: str) -> None:
    from tools.mcp_tool_lifecycle import shutdown_mcp_servers
    shutdown_mcp_servers(names={name})


def load_plugins() -> None:
    from hermes_cli.plugins import discover_plugins
    discover_plugins()
