"""The JSON lines the runtime server and an agent process exchange (`AgentCommand` and `AgentProcessEvent` in
src/types.ts)."""
from collections.abc import Callable
from typing import Any, Literal, NotRequired, TypedDict

JsonObject = dict[str, Any]


class Image(TypedDict):
    mediaType: str
    data: str


class ClientTool(TypedDict):
    name: str
    description: str
    inputSchema: NotRequired[JsonObject]


class McpServer(TypedDict):
    url: str
    headers: dict[str, str]


class Limits(TypedDict):
    maxTurns: int
    maxBudgetUsd: float


class Settings(TypedDict):
    homes: str
    skillsDir: str
    compressionThresholdTokens: int
    warm: bool


class InitCommand(TypedDict):
    type: Literal['init']
    credentials: dict[str, str]
    settings: Settings


class TurnCommand(TypedDict):
    type: Literal['turn']
    text: str
    images: NotRequired[list[Image]]
    systemPrompt: str
    effort: NotRequired[str]
    model: NotRequired[str]
    cwd: str
    userDir: NotRequired[str]
    outputSchema: NotRequired[JsonObject]
    clientTools: NotRequired[list[ClientTool]]
    datagrokMcp: NotRequired[McpServer]
    limits: Limits
    restore: NotRequired[str]


class ContinueCommand(TypedDict):
    type: Literal['continue']
    text: str


class CancelCommand(TypedDict):
    type: Literal['cancel']


class BrowserResultCommand(TypedDict):
    type: Literal['browser_result']
    id: str
    value: NotRequired[Any]
    error: NotRequired[str]


Command = TurnCommand | ContinueCommand | CancelCommand | BrowserResultCommand


class Metrics(TypedDict):
    inputTokens: int
    outputTokens: int
    cacheReadTokens: int
    cacheCreationTokens: int
    costUsd: float
    numTurns: int
    durationMs: int
    durationApiMs: int


class ReadyEvent(TypedDict):
    type: Literal['ready']


class FatalEvent(TypedDict):
    type: Literal['fatal']
    message: str


class SessionEvent(TypedDict):
    type: Literal['session']
    id: str


class TextEvent(TypedDict):
    type: Literal['text']
    text: str


class ActivityEvent(TypedDict):
    type: Literal['activity']


class ToolStartEvent(TypedDict):
    type: Literal['tool_start']
    name: str
    args: JsonObject


class ToolEndEvent(TypedDict):
    type: Literal['tool_end']
    name: str
    args: JsonObject
    result: str


class BrowserCallEvent(TypedDict):
    type: Literal['browser_call']
    id: str
    name: str
    input: JsonObject
    timeoutMs: int


class StoppedEvent(TypedDict):
    type: Literal['stopped']
    text: str
    structuredOutput: NotRequired[Any]
    metrics: Metrics


class LimitEvent(TypedDict):
    type: Literal['limit']
    kind: Literal['turns', 'budget']


class ErrorEvent(TypedDict):
    type: Literal['error']
    message: str


class CancelledEvent(TypedDict):
    type: Literal['cancelled']


class SessionLostEvent(TypedDict):
    type: Literal['session_lost']


Event = (ReadyEvent | FatalEvent | SessionEvent | TextEvent | ActivityEvent | ToolStartEvent | ToolEndEvent |
         BrowserCallEvent | StoppedEvent | LimitEvent | ErrorEvent | CancelledEvent | SessionLostEvent)
Emit = Callable[[Event], None]
