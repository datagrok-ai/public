"""Tools that run in the user's Datagrok tab. Each call goes to Node, which forwards it to the browser and returns its
answer. They are native Hermes tools, so their results reach the model as they are."""
import json
import re
from collections.abc import Callable
from typing import Any

from jsonschema import Draft7Validator

from .hermes import ToolHandler, ToolRegistry, ToolSpec
from .protocol import ClientTool, JsonObject

TOOLSET = 'grokky-browser'
VIEW_TOOLSET = 'grokky-view'

# Seconds. An unanswered question must not hold the turn (and the workspace sync behind it) open.
EXEC_TIMEOUT, VERIFY_TIMEOUT, ENTITIES_TIMEOUT, VIEW_TIMEOUT, ASK_TIMEOUT = 300, 60, 30, 120, 60

ENTITY_TYPES = {'dataquery': 'query', 'dataconnection': 'connection', 'usergroup': 'group', 'fileinfo': 'file'}

CallBrowser = Callable[[str, JsonObject, int], Any]


def _string(description: str) -> JsonObject:
    return {'type': 'string', 'description': description}


EXEC: ToolSpec = {
    'name': 'datagrok_exec',
    'description': (
        'Run JavaScript in the Datagrok tab to perform an action (add viewer, filter, open file, upload data, …). '
        'Returns {success, returnValue?, verified?, error?}. ALWAYS provide `verify` with state-changing code — it '
        'runs right after your code in the same round-trip and its {passed, observed} comes back as `verified`, so '
        'no separate datagrok_verify call is needed when it passes. For informational questions ("how do I…", '
        '"what is…") answer in plain text — do NOT call this tool.'),
    'parameters': {'type': 'object', 'required': ['code'], 'properties': {
        'code': _string('Async JS (await works). Globals: grok, ui, DG, view, t. Return a plain object confirming the '
                        'action (shape per the datagrok-exec skill), or an HTMLElement only to render output in chat.'),
        'verify': {'type': 'object', 'required': ['assertion', 'description'],
                   'description': 'Required for any state-changing code; omit only when the code is a pure '
                                  'read/render.',
                   'properties': {
                       'assertion': _string(
                           'Async JS that re-reads the state your code changed and MUST `return` the values it '
                           'read (truthy = verified; a bare `true` is rejected). Runs in a fresh scope (globals grok, '
                           "ui, DG, view, t) — it cannot see your code's variables, so re-derive everything from live "
                           'state.'),
                       'description': _string('Short human-readable statement of what is being verified.')}}}},
}

VERIFY: ToolSpec = {
    'name': 'datagrok_verify',
    'description': (
        'Confirm a previous datagrok_exec / MCP action actually took effect, by re-reading LIVE state. Runs in a fresh '
        'scope (globals: grok, ui, DG, view, t) — it CANNOT see variables from your action code, so it must re-derive '
        'from t/view/grok. Returns {passed, observed?, error?}. After any action you MUST call this before reporting '
        'success; passed:false means it did NOT work — fix it, do not report done.'),
    'parameters': {'type': 'object', 'required': ['assertion', 'description'], 'properties': {
        'assertion': _string('Async JS (await works) that MUST `return` the values it read; truthy = verified. '
                             'Re-read whatever you changed from t/view/grok. A bare `true` or any literal is '
                             'rejected.'),
        'description': _string('Short human-readable statement of what is being verified.')}},
}

SHOW_ENTITIES: ToolSpec = {
    'name': 'datagrok_show_entities',
    'description': (
        'Display Datagrok entities (files, scripts, queries, connections, projects, spaces, groups, users) as '
        'interactive cards in the chat. ALWAYS call this tool when you have entity data from an MCP result — never '
        'list entities as plain text, markdown links, or a fenced block.'),
    'parameters': {'type': 'object', 'required': ['entities'], 'properties': {'entities': {
        'type': 'array',
        'description': (
            'Array of entity refs. Each entry: type (file|script|query|connection|project|space|group|user) or the raw '
            '#type string from the MCP result — the tool normalizes it; name (for users, use friendlyName); id '
            '(required for all types except file); connector + path (required for files; connector is the connection '
            'name, e.g. "System:DemoFiles"); isDirectory (optional, files only). Pass all values from the MCP result '
            'exactly — never invent.'),
        'items': {'type': 'object', 'required': ['name'], 'properties': {
            '#type': {'type': 'string'}, 'type': {'type': 'string'}, 'name': {'type': 'string'},
            'id': {'type': 'string'}, 'connector': {'type': 'string'}, 'path': {'type': 'string'},
            'isDirectory': {'type': 'boolean'}}}}}},
}

ASK_USER: ToolSpec = {
    'name': 'AskUserQuestion',
    'description': ('Ask the user a multiple-choice question when the intent of their request is genuinely unclear. '
                    "Returns the user's answers."),
    'parameters': {'type': 'object', 'required': ['questions'], 'properties': {'questions': {
        'type': 'array', 'description': '1-4 questions to ask.',
        'items': {'type': 'object', 'required': ['question', 'options'], 'properties': {
            'question': _string('The complete question to ask the user.'),
            'header': _string('A very short label for the question (at most 12 characters).'),
            'options': {'type': 'array', 'description': '2-4 distinct choices.', 'items': {
                'type': 'object', 'required': ['label'], 'properties': {
                    'label': _string('The choice as shown to the user (1-5 words).'),
                    'description': _string('What this choice means.')}}},
            'multiSelect': {'type': 'boolean', 'description': 'Allow the user to pick more than one option.'}}}}}},
}

# What the model API accepts as a tool name.
TOOL_NAME = re.compile(r'[A-Za-z0-9_-]{1,64}')


def _entity(e: JsonObject) -> JsonObject:
    raw = str(e.get('#type') or e.get('type') or '').lower()
    out: JsonObject = {'type': ENTITY_TYPES.get(raw, raw), 'name': e.get('name')}
    out.update({k: e[k] for k in ('id', 'connector', 'path') if e.get(k)})
    if 'isDirectory' in e:
        out['isDirectory'] = e['isDirectory']
    return out


def _entities(args: JsonObject) -> JsonObject:
    return {'entities': [_entity(e) for e in args.get('entities') or []]}


def _invalid(spec: ToolSpec, args: JsonObject) -> str | None:
    errors = list(Draft7Validator(spec['parameters']).iter_errors(args))
    return f"Input validation error: Invalid arguments for tool {spec['name']}: {errors[0].message}" if errors else None


def _relay(spec: ToolSpec, call: CallBrowser, timeout: int, failure_key: str,
           prepare: Callable[[JsonObject], JsonObject] = lambda args: args) -> ToolHandler:
    def handler(args: JsonObject, **_: Any) -> str:
        try:
            invalid = _invalid(spec, args)
            if invalid:
                return json.dumps({failure_key: False, 'error': invalid})
            return json.dumps(call(spec['name'], prepare(args), timeout))
        except Exception as e:
            return json.dumps({failure_key: False, 'error': str(e)})
    return handler


def _ask_user(call: CallBrowser) -> ToolHandler:
    def handler(args: JsonObject, **_: Any) -> str:
        invalid = _invalid(ASK_USER, args)
        if invalid:
            return invalid
        try:
            answers = (call(ASK_USER['name'], args, ASK_TIMEOUT) or {}).get('answers') or {}
        except Exception:
            return 'The user did not answer the question in time. Proceed without the answer or end the turn.'
        given = ', '.join(f'"{q}"="{a}"' for q, a in answers.items())
        return f'Your questions have been answered: {given}. You can now continue with these answers in mind.'
    return handler


def register(registry: ToolRegistry, call: CallBrowser) -> None:
    """The browser tools every session has."""
    tools = [(EXEC, _relay(EXEC, call, EXEC_TIMEOUT, 'success')),
             (VERIFY, _relay(VERIFY, call, VERIFY_TIMEOUT, 'passed')),
             (SHOW_ENTITIES, _relay(SHOW_ENTITIES, call, ENTITIES_TIMEOUT, 'success', _entities)),
             (ASK_USER, _ask_user(call))]
    for spec, handler in tools:
        registry.register(name=spec['name'], toolset=TOOLSET, schema=spec, handler=handler)


def view_tool_specs(client_tools: list[ClientTool] | None) -> list[ToolSpec]:
    """The browser's view tools that can be offered: well-named and described."""
    specs: list[ToolSpec] = []
    for t in client_tools or []:
        name = t.get('name') if isinstance(t, dict) else None
        if isinstance(name, str) and TOOL_NAME.fullmatch(name) and t.get('description'):
            specs.append({'name': name, 'description': t['description'],
                          'parameters': {'type': 'object', 'properties': {}, **(t.get('inputSchema') or {})}})
    return specs


def register_view_tools(registry: ToolRegistry, call: CallBrowser, specs: list[ToolSpec], previous: list[str]) \
        -> list[str]:
    """Replaces the previous turn's view tools with these; returns the names it registered. A name another tool has
    is skipped: Hermes refuses it, and deregistering it next turn would remove that tool."""
    for name in previous:
        registry.deregister(name)
    registered = []
    for spec in specs:
        if registry.get_entry(spec['name']) is None:
            registry.register(name=spec['name'], toolset=VIEW_TOOLSET, schema=spec,
                              handler=_relay(spec, call, VIEW_TIMEOUT, 'success'))
            registered.append(spec['name'])
    return registered
