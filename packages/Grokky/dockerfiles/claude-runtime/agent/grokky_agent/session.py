"""One Grokky chat on Hermes' AIAgent: the turns Node sends, and the Hermes hooks of this process."""
import json
import re
import threading
import time
import uuid
from dataclasses import dataclass, field
from typing import Any, Literal

from . import browser_tools, guard, hermes, models, structured
from .models import Usage
from .protocol import BrowserResultCommand, ContinueCommand, Emit, JsonObject, McpServer, Metrics, TurnCommand
from .providers import Provider

MAX_STRUCTURED_ATTEMPTS = 5
HEARTBEAT_SECONDS = 30
# Node answers every browser call within its own timeout; this only frees a thread whose call Node dropped.
BROWSER_GRACE_SECONDS = 30
# Hermes stamps its activity clock with this for every stream event (`_count_chunk`, agent/chat_completion_helpers.py).
STREAM_EVENT = 'receiving stream response'
# Bounds the connect, the MCP handshake and the tool listing: Node's idle watchdog ends a silent turn after 90 s.
MCP_CONNECT_TIMEOUT = 30
MCP_SERVER = 'datagrok'
# Hermes reports a failed call as "Error code: 400 - {...'message': '<provider text>'}" (a text with an apostrophe in
# double quotes, JSON in double quotes), or with the HTML page a proxy or gateway answered instead; the user sees the
# text, or the page's title. The title is trimmed in code: trimming it in the pattern backtracks on an unclosed one.
PROVIDER_MESSAGE = re.compile(r"""['"]message['"]: (?:'((?:[^'\\]|\\.)*)'|"((?:[^"\\]|\\.)*)")|"""
                              r'(?i:<title>([^<]*)</title>)')
# api_request_error's reason when the provider's safety filter stops an answer; its message is the cut-off answer.
SAFETY_FILTER = 'content_policy_blocked'
SAFETY_FILTER_MESSAGE = "The model provider's safety filter stopped this answer. Try rewording the request."
BASE_TOOLSETS = ['file', 'terminal', 'skills', 'web', browser_tools.TOOLSET]
RETRY_STRUCTURED = ('Your reply did not match the required JSON schema: {}. '
                    'Reply again with only the corrected JSON object.')

current: 'Session | None' = None


@dataclass
class Tokens:
    input: int = 0
    output: int = 0
    cache_read: int = 0
    cache_write: int = 0

    def add(self, usage: Usage) -> None:
        self.input += usage.get('input_tokens') or 0
        self.output += usage.get('output_tokens') or 0
        self.cache_read += usage.get('cache_read_tokens') or 0
        self.cache_write += usage.get('cache_write_tokens') or 0


class Turn:
    """One browser message, with its gate follow-ups and structured-output retries."""

    def __init__(self, cmd: TurnCommand) -> None:
        self.system_prompt = cmd['systemPrompt']
        self.schema = cmd.get('outputSchema')
        self.cwd = cmd['cwd']
        self.user_dir = cmd.get('userDir')
        self.max_calls = cmd['limits']['maxTurns']
        self.max_cost = cmd['limits']['maxBudgetUsd']
        self.started = time.time()
        self.api_calls = 0
        self.api_seconds = 0.0
        self.cost = 0.0
        self.tokens = Tokens()
        self.limit: Literal['turns', 'budget'] | None = None
        self.api_error: str | None = None
        self.cancelled = False
        self.segment_start = True
        self.in_flight: set[str] = set()

    def metrics(self) -> Metrics:
        return {'inputTokens': self.tokens.input, 'outputTokens': self.tokens.output,
                'cacheReadTokens': self.tokens.cache_read, 'cacheCreationTokens': self.tokens.cache_write,
                'costUsd': self.cost, 'numTurns': self.api_calls,
                'durationMs': round((time.time() - self.started) * 1000),
                'durationApiMs': round(self.api_seconds * 1000)}


@dataclass
class BrowserCall:
    answered: threading.Event = field(default_factory=threading.Event)
    value: Any = None
    error: str | None = None

    def answer(self, value: Any = None, error: str | None = None) -> None:
        self.value, self.error = value, error
        self.answered.set()


class Heartbeat:
    """Tells Node a turn is busy, at most every 30 s, while a tool or auxiliary call runs or the model streams; Hermes'
    other stamps on its activity clock are silence and must still trip Node's idle watchdog."""

    def __init__(self, emit: Emit) -> None:
        self._emit = emit
        self._last_at: float | None = None
        self._busy = False
        self._sent = float('-inf')

    @staticmethod
    def streamed(last_at: float | None, activity: hermes.Activity) -> bool:
        return activity.get('last_activity_at') != last_at and activity.get('last_activity_description') == STREAM_EVENT

    def beat(self, now: float, turn: Turn | None, agent: hermes.Agent | None) -> None:
        if turn is None or agent is None:
            self._busy = False
            return
        activity = agent.get_activity_summary()
        self._busy = self._busy or bool(turn.in_flight) or self.streamed(self._last_at, activity)
        self._last_at = activity.get('last_activity_at')
        if self._busy and now - self._sent >= HEARTBEAT_SECONDS:
            self._emit({'type': 'activity'})
            self._busy, self._sent = False, now


class Session:
    def __init__(self, emit: Emit, provider: Provider, homes: str, store: hermes.SessionStore,
                 registry: hermes.ToolRegistry, spillover: str) -> None:
        self.emit = emit
        self.provider = provider
        self.homes = homes
        self.store = store
        self.registry = registry
        self.spillover = spillover
        self.key = uuid.uuid4().hex
        self.agent: hermes.Agent | None = None
        self.model = provider.default
        self.history: list[JsonObject] = []
        self.hermes_session: str | None = None
        # The turn whose run is in progress (the hooks act for it), and the latest one received, which gate feedback
        # continues and a cancel stops even before it runs.
        self.turn: Turn | None = None
        self.latest: Turn | None = None
        self.view_tools: list[str] = []
        self.view_specs: list[hermes.ToolSpec] = []
        self.mcp: McpServer | None = None
        # The stored session this process reloads, until a turn has restored it.
        self.restore: str | None = None
        self.calls: dict[str, BrowserCall] = {}
        self.hooks_registered = False
        browser_tools.register(registry, self.call_browser)
        threading.Thread(target=self._beat_every_second, daemon=True).start()

    def register_hooks(self, ctx: hermes.PluginContext) -> None:
        hooks: dict[str, hermes.Hook] = {
            'pre_tool_call': self.pre_tool_call, 'post_tool_call': self.post_tool_call,
            'pre_api_request': self.pre_api_request, 'post_api_request': self.post_api_request,
            'api_request_error': self.api_request_error, 'pre_auxiliary_call': self.pre_auxiliary_call,
            'post_auxiliary_call': self.post_auxiliary_call}
        for name, callback in hooks.items():
            ctx.register_hook(name, callback)
        ctx.register_middleware('llm_request', self.llm_request)
        self.hooks_registered = True

    def pre_tool_call(self, tool_name: str, args: JsonObject | None = None, tool_call_id: str = '',
                      **_: Any) -> hermes.ToolVerdict | None:
        t = self.turn
        if t is None or t.cancelled:
            return {'action': 'block', 'message': 'No turn is in progress.'}
        args = args or {}
        self.emit({'type': 'tool_start', 'name': tool_name, 'args': args})
        why = guard.violation(tool_name, args, t.user_dir, self.homes, self.spillover)
        if why:
            return {'action': 'block', 'message': why}
        t.in_flight.add(tool_call_id or tool_name)
        path = guard.absolute_edit_path(tool_name, args, t.cwd)
        return {'action': 'modify', 'args': {'path': path}} if path else None

    def post_tool_call(self, tool_name: str, args: JsonObject | None = None, result: Any = None,
                       status: str | None = None, tool_call_id: str = '', **_: Any) -> None:
        t = self.turn
        if t is None:
            return
        t.in_flight.discard(tool_call_id or tool_name)
        t.segment_start = True
        # Hermes reports a call the guard refused too; it never ran.
        if status != 'blocked':
            self.emit({'type': 'tool_end', 'name': tool_name, 'args': args or {},
                       'result': result if isinstance(result, str) else json.dumps(result, default=str)})

    def pre_api_request(self, **_: Any) -> None:
        # Every attempt, retries included: the error of the one before no longer stands.
        if self.turn:
            self.turn.api_error = None

    def post_api_request(self, usage: Usage | None = None, api_duration: float | None = None,
                         finish_reason: str | None = None, **_: Any) -> None:
        t = self.turn
        if t is None:
            return
        t.api_calls += 1
        self._account(t, usage or {}, api_duration)
        # The reply's tools must not run when no call is left to answer them with.
        if t.api_calls >= t.max_calls and finish_reason in ('tool_calls', 'length'):
            self._stop(t, 'turns')

    def api_request_error(self, status_code: int | None = None, error: hermes.ApiError | None = None,
                          reason: str | None = None, **_: Any) -> None:
        t = self.turn
        if t is None:
            return
        if reason == SAFETY_FILTER:
            t.api_error = SAFETY_FILTER_MESSAGE
            return
        raw = (error and error['message']) or reason or 'The model request failed.'
        m = PROVIDER_MESSAGE.search(raw)
        message = (m.group(1) or m.group(2) or m.group(3) or '').strip() if m else ''
        text = f'API Error: {status_code} {message}' if message and status_code else raw
        t.api_error = text + self.provider.error_hint(status_code)

    def pre_auxiliary_call(self, api_request_id: str = '', **_: Any) -> None:
        if self.turn:
            self.turn.in_flight.add(f'aux:{api_request_id}')

    def post_auxiliary_call(self, usage: Usage | None = None, api_duration: float | None = None,
                            api_request_id: str = '', **_: Any) -> None:
        t = self.turn
        if t is None:
            return
        t.in_flight.discard(f'aux:{api_request_id}')
        self._account(t, usage or {}, api_duration)

    def llm_request(self, request: JsonObject | None = None, api_mode: str = '', **_: Any) -> hermes.RequestOverride:
        t = self.turn
        if t is None or request is None:
            return {'request': self.provider.refused_request(api_mode), 'source': 'grokky'}
        # Hermes also calls the model on its own: after an empty or thinking-only reply, for a summary at its own limit.
        if t.api_calls >= t.max_calls:
            self._stop(t, 'turns')
            return {'request': self.provider.refused_request(api_mode), 'source': 'grokky'}
        # Hermes sends its own request when this raises, so a patch that fails sends the refusal.
        try:
            patched = self.provider.patch_request(api_mode, request, t.system_prompt, t.schema)
        except Exception:
            patched = self.provider.refused_request(api_mode)
        return {'request': patched, 'source': 'grokky'}

    def receive(self, cmd: TurnCommand | ContinueCommand) -> Turn:
        """The turn a `turn` or `continue` command runs, taken when the command arrives."""
        if cmd['type'] == 'turn':
            self.latest = Turn(cmd)
        assert self.latest is not None, 'a continue follows its turn'
        return self.latest

    def start_turn(self, cmd: TurnCommand, t: Turn) -> None:
        self.restore = cmd.get('restore') or self.restore
        model = cmd.get('model') or self.provider.default
        if model not in self.provider.models:
            self.emit({'type': 'error',
                       'message': f'Unknown model "{model}": use one of {", ".join(self.provider.models)}'})
            return
        if cmd.get('images') and not self.provider.images(model):
            self.emit({'type': 'error', 'message': f'The model {model} does not accept images'})
            return
        problem = structured.schema_error(t.schema) if t.schema is not None else None
        if problem:
            self.emit({'type': 'error', 'message': f'The output schema cannot be used: {problem}'})
            return
        agent = self._prepare(cmd, model)
        if agent is None:
            self.emit({'type': 'session_lost'})
            return
        effort = cmd.get('effort')
        agent.reasoning_config = {'enabled': True, 'effort': effort} if effort else {'enabled': False}
        self._run(agent, t, self._user_content(cmd))

    def continue_turn(self, cmd: ContinueCommand, t: Turn) -> None:
        assert self.agent is not None, 'a continue follows its turn'
        self._run(self.agent, t, cmd['text'])

    def cancel(self) -> None:
        t = self.latest
        if t is None or t.cancelled:
            return
        t.cancelled = True
        if self.turn is t and self.agent:
            self.agent.interrupt('cancelled by the user', hard_cancel=True)
        for call in list(self.calls.values()):
            call.answer(error='aborted')

    def browser_answer(self, cmd: BrowserResultCommand) -> None:
        call = self.calls.get(cmd['id'])
        if call:
            call.answer(cmd.get('value'), cmd.get('error'))

    def call_browser(self, name: str, args: JsonObject, timeout_s: int) -> Any:
        call_id = uuid.uuid4().hex
        call = self.calls[call_id] = BrowserCall()
        self.emit({'type': 'browser_call', 'id': call_id, 'name': name, 'input': args, 'timeoutMs': timeout_s * 1000})
        answered = call.answered.wait(timeout_s + BROWSER_GRACE_SECONDS)
        del self.calls[call_id]
        if not answered:
            raise RuntimeError('The browser did not answer in time.')
        if call.error:
            raise RuntimeError(call.error)
        return call.value

    def build_agent(self, model: str, hermes_session: str | None = None) -> hermes.Agent:
        self.hermes_session = hermes_session or uuid.uuid4().hex
        self.agent = hermes.new_agent(model, self.provider.connection(model), self.hermes_session, self.store,
                                      self._toolsets(), self._on_delta)
        self.model = model
        self.emit({'type': 'session', 'id': self.hermes_session})
        return self.agent

    def _prepare(self, cmd: TurnCommand, model: str) -> hermes.Agent | None:
        """The agent for this turn's model and stored session, with this turn's tools; None when the stored session
        to continue is gone."""
        if self.restore:
            agent = self.build_agent(model, self.restore)
            self.history = self.store.get_messages_as_conversation(self.restore, repair_alternation=True)
            if not self.history:
                return None
            self.restore = None
        elif self.agent is None or model != self.model:
            agent = self.build_agent(model, self.hermes_session)
        else:
            agent = self.agent
        hermes.set_terminal_cwd(self.key, cmd['cwd'])
        changed = self._set_view_tools(browser_tools.view_tool_specs(cmd.get('clientTools')))
        changed = self._set_mcp(cmd.get('datagrokMcp')) or changed
        if changed:
            agent.enabled_toolsets = self._toolsets()
            agent.tools = hermes.tool_definitions(agent.enabled_toolsets)
            agent.valid_tool_names = {t['function']['name'] for t in agent.tools}
        return agent

    def _toolsets(self) -> list[str]:
        return BASE_TOOLSETS + ([browser_tools.VIEW_TOOLSET] if self.view_tools else []) + \
            ([f'mcp-{MCP_SERVER}'] if self.mcp else [])

    def _set_view_tools(self, specs: list[hermes.ToolSpec]) -> bool:
        if specs == self.view_specs:
            return False
        self.view_tools = browser_tools.register_view_tools(self.registry, self.call_browser, specs, self.view_tools)
        self.view_specs = specs
        return True

    def _set_mcp(self, server: McpServer | None) -> bool:
        """The Datagrok MCP server, reached with this user's key: Hermes keeps one connection per server name, so a
        changed key or URL replaces the connection."""
        if server == self.mcp:
            return False
        if self.mcp:
            hermes.disconnect_mcp(MCP_SERVER)
        self.mcp = None
        if server:
            names = hermes.connect_mcp(MCP_SERVER, {'url': server['url'], 'headers': server['headers'],
                                                    'connect_timeout': MCP_CONNECT_TIMEOUT})
            if any(n.startswith(f'mcp__{MCP_SERVER}__') for n in names):
                self.mcp = server
            else:
                # Hermes parks a server it could not reach and reconnects it only in the background when asked again,
                # so the next turn would still start without its tools: dropped, it is connected afresh.
                hermes.disconnect_mcp(MCP_SERVER)
        return True

    def _run(self, agent: hermes.Agent, t: Turn, content: str | list[JsonObject]) -> None:
        # An interrupt that landed after the previous run ended is stale (Hermes clears its own in the finalizer).
        agent.clear_interrupt()
        self.turn = t
        try:
            if t.cancelled:
                self.emit({'type': 'cancelled'})
            else:
                self._prompt(agent, t, content)
        finally:
            self.turn = None

    def _prompt(self, agent: hermes.Agent, t: Turn, content: str | list[JsonObject]) -> None:
        problem = ''
        for _ in range(MAX_STRUCTURED_ATTEMPTS):
            t.segment_start = True
            # One above what is left, so the request patch, not Hermes' own limit, ends the turn at the limit.
            agent.max_iterations = max(1, t.max_calls - t.api_calls) + 1
            result = agent.run_conversation(user_message=content, conversation_history=self.history, task_id=self.key)
            if isinstance(result.get('messages'), list):
                self.history = result['messages']
            if agent.session_id != self.hermes_session:
                self.hermes_session = agent.session_id
                self.emit({'type': 'session', 'id': agent.session_id})
            if t.cancelled:
                self.emit({'type': 'cancelled'})
                return
            if t.limit or str(result.get('turn_exit_reason') or '').startswith('max_iterations'):
                self.emit({'type': 'limit', 'kind': t.limit or 'turns'})
                return
            if t.api_error or result.get('failed'):
                self.emit({'type': 'error',
                           'message': t.api_error or result.get('error') or 'The model request failed.'})
                return
            text = result.get('final_response') or ''
            if t.schema is None:
                if text:
                    self.emit({'type': 'stopped', 'text': text, 'metrics': t.metrics()})
                else:
                    self.emit({'type': 'error', 'message': 'The assistant stopped without an answer.'})
                return
            value, error = structured.check(t.schema, text)
            if error is None:
                self.emit({'type': 'stopped', 'text': text, 'structuredOutput': value, 'metrics': t.metrics()})
                return
            problem = error
            content = RETRY_STRUCTURED.format(problem)
        self.emit({'type': 'error', 'message': 'The model did not return valid structured output after '
                                               f'{MAX_STRUCTURED_ATTEMPTS} attempts: {problem}'})

    def _account(self, t: Turn, usage: Usage, api_duration: float | None) -> None:
        t.tokens.add(usage)
        t.api_seconds += api_duration or 0.0
        t.cost += models.cost(self.provider.models[self.model], usage)
        if t.cost > t.max_cost:
            self._stop(t, 'budget')

    def _stop(self, t: Turn, limit: Literal['turns', 'budget']) -> None:
        if t.limit is None and self.agent:
            t.limit = limit
            self.agent.interrupt(f'{limit} limit', hard_cancel=True)

    def _on_delta(self, text: str | None) -> None:
        t = self.turn
        if t is None or not text:
            return
        # Hermes opens the first text after a tool round with a paragraph break for its own display.
        if t.segment_start:
            text = text.lstrip('\n')
            if not text:
                return
            t.segment_start = False
        # A structured reply is checked whole, and a failed call's prose is not the answer.
        if t.schema is not None or t.api_error:
            self.emit({'type': 'activity'})
        else:
            self.emit({'type': 'text', 'text': text})

    @staticmethod
    def _user_content(cmd: TurnCommand) -> str | list[JsonObject]:
        images = cmd.get('images')
        if not images:
            return cmd['text']
        return [{'type': 'image_url', 'image_url': {'url': f"data:{i['mediaType']};base64,{i['data']}"}}
                for i in images] + [{'type': 'text', 'text': cmd['text']}]

    def _beat_every_second(self) -> None:
        heartbeat = Heartbeat(self.emit)
        # Hermes stamps its own waits on the clock the stream events use, every 30 s: sampling must be finer.
        while True:
            time.sleep(1)
            heartbeat.beat(time.monotonic(), self.turn, self.agent)
