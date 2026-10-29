"""Tool calls a model writes as its answer instead of making them (Llama on most hosts): `name(k=v, ...)`,
`name(op:'x', args:{...})`, `{"name": ..., "parameters": {...}}` or `<function=name>{...}</function>`. A reply that is
one such call to a tool of the session is run as that call; anything else is an answer."""
import ast
import re
from typing import Any

CALL = re.compile(r'^(?:mcp__[a-z0-9]+__)?([A-Za-z_][A-Za-z0-9_]*)\s*\((.*)\)\s*;?$', re.S)
TAGGED = re.compile(r'^<function\s*(?:=|name\s*=\s*)"?([A-Za-z_][A-Za-z0-9_]*)"?\s*>(.*?)(?:</function>)?$', re.S)
BARE_KEY = re.compile(r'([{,]\s*)([A-Za-z_][A-Za-z0-9_]*)\s*:')
JSON_LITERALS = {'true': True, 'false': False, 'null': None}


def _literal(node: ast.expr) -> Any:
    if isinstance(node, ast.Name) and node.id in JSON_LITERALS:
        return JSON_LITERALS[node.id]
    if isinstance(node, ast.Dict):
        return {_literal(k): _literal(v) for k, v in zip(node.keys, node.values, strict=True) if k is not None}
    if isinstance(node, ast.List):
        return [_literal(e) for e in node.elts]
    return ast.literal_eval(node)


def _object(text: str) -> dict[str, Any] | None:
    """A JSON object, or a JS/Python-style one with bare keys and single quotes."""
    for candidate in (text, BARE_KEY.sub(r'\1"\2":', text)):
        try:
            value = _literal(ast.parse(candidate, mode='eval').body)
        except (SyntaxError, ValueError):
            continue
        if isinstance(value, dict):
            return value
    return None


def _arguments(body: str) -> dict[str, Any] | None:
    body = body.strip()
    if not body:
        return {}
    try:
        call = ast.parse(f'f({body})', mode='eval').body
        if isinstance(call, ast.Call) and not call.args and all(k.arg for k in call.keywords):
            return {str(k.arg): _literal(k.value) for k in call.keywords}
    except (SyntaxError, ValueError):
        pass
    return _object('{' + body + '}')


def parse(text: str, tools: list[str]) -> tuple[str, dict[str, Any]] | None:
    """(tool, arguments) when the reply is one call to one of `tools`, alone or ending a few words of prose
    ("Let me try again. name(...)"), else None."""
    text = re.sub(r'^```\w*\s*|\s*```$', '', text.strip())
    if (call := _parse_call(text, tools)) is not None:
        return call
    names = '|'.join(re.escape(re.sub(r'^mcp__.+?__', '', t)) for t in tools)
    starts = [m.start() for m in re.finditer(rf'^\s*(?:mcp__[a-z0-9]+__)?(?:{names})\s*\(', text, re.M)]
    return _parse_call(text[starts[-1]:].strip(), tools) if starts else None


def _parse_call(text: str, tools: list[str]) -> tuple[str, dict[str, Any]] | None:
    name, args = '', None
    if m := CALL.match(text):
        name, args = m.group(1), _arguments(m.group(2))
    elif m := TAGGED.match(text):
        name, args = m.group(1), _object(m.group(2).strip() or '{}')
    elif text.startswith('{') and (obj := _object(text)) is not None:
        name = str(obj.get('name') or '')
        params = obj.get('parameters', obj.get('arguments', {}))
        args = params if isinstance(params, dict) else None
    if not name or args is None:
        return None
    match = next((t for t in tools if t == name), None) or next((t for t in tools if t.endswith('_' + name)), None)
    return (match, args) if match else None


def result_message(tool: str, result: str) -> str:
    return (f'[The {tool} call you wrote ran; its result follows. Next time make the call through the tool interface '
            f'instead of writing it.]\n{result}')
