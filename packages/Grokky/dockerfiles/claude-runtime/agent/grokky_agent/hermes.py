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
    # On the Anthropic Messages API Hermes sends a string key as x-api-key (send_keys_by_type) and a callable's as
    # Bearer, except to api.anthropic.com, which gets a Console key as x-api-key either way.
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


class Usage(TypedDict, total=False):
    input_tokens: int | None
    output_tokens: int | None
    cache_read_tokens: int | None
    cache_write_tokens: int | None


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


def send_keys_by_type() -> None:
    """Hermes would send a string key as Bearer to the hosts on its list, Azure's among them, which documents keys in
    x-api-key or api-key and Bearer for Entra ID: without the list the key's type decides (Connection), on every
    Anthropic client of the process."""
    import agent.anthropic_adapter as anthropic_adapter
    anthropic_adapter._requires_bearer_auth = lambda base_url: False


def use_bedrock_key(key: str) -> None:
    """Hermes reaches Bedrock's Converse API through boto3, which takes a Bedrock API key only from the process
    environment, where the session's tools could read it from /proc. The clients are built here instead: unsigned,
    with the key added as a Bearer header, so boto3 never walks its credential chain either."""
    import agent.bedrock_adapter as bedrock_adapter
    import boto3
    import botocore
    from botocore.config import Config

    clients: dict[tuple[str, str], Any] = {}

    def add_key(request: Any, **_: Any) -> None:
        request.headers['Authorization'] = f'Bearer {key}'

    def client(service: str, region: str) -> Any:
        if (service, region) not in clients:
            # boto3 stops reading after 60 s, under Kimi's long replies; Hermes retries and spots stale streams.
            c = boto3.client(service, region_name=region, config=Config(
                signature_version=botocore.UNSIGNED, read_timeout=600, retries={'max_attempts': 1}))
            c.meta.events.register(f'before-send.{service}.*', add_key)
            clients[service, region] = c
        return clients[service, region]

    bedrock_adapter._get_bedrock_runtime_client = lambda region: client('bedrock-runtime', region)
    bedrock_adapter._get_bedrock_control_client = lambda region: client('bedrock', region)
    bedrock_adapter.resolve_bedrock_bearer_token = lambda env=None: key
    # Hermes learns a window by sending ~1.3M tokens to provoke a length error: over 30 s for some models (Kimi,
    # Qwen), and a billed call for any that accepts it. Its table, or 128k, serves instead.
    bedrock_adapter.probe_bedrock_context_length = lambda model_id, region: None


def usage_cost(model: str, usage: Usage, billing: str, base_url: str) -> float | None:
    """USD for one call from Hermes' price list for the provider `billing` names, or None when it has no price."""
    from agent.usage_pricing import CanonicalUsage, estimate_usage_cost
    canonical = CanonicalUsage(input_tokens=usage.get('input_tokens') or 0,
                               output_tokens=usage.get('output_tokens') or 0,
                               cache_read_tokens=usage.get('cache_read_tokens') or 0,
                               cache_write_tokens=usage.get('cache_write_tokens') or 0)
    amount = estimate_usage_cost(model, canonical, provider=billing, base_url=base_url).amount_usd
    return None if amount is None else float(amount)


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


def describe_skills(skills: list[tuple[str, str]]) -> None:
    """Puts the index of Grokky's skills into the skill_view tool's own description: models that weigh tool
    definitions over a long system prompt (Qwen) otherwise act without loading the skill the request matches."""
    from tools.skills_tool import SKILL_VIEW_SCHEMA
    SKILL_VIEW_SCHEMA['description'] = (
        'Load the Datagrok skill whose description matches the request BEFORE acting: it gives the exact API calls '
        'and conventions for the task. Skills:\n' + '\n'.join(f'- {name}: {desc}' for name, desc in skills))



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
