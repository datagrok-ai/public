"""The request formats a model is called in, by the api_mode Hermes names them with: where a turn's system prompt and
output schema go on a request, and the request sent instead of one the session cannot vouch for."""
import copy
from collections.abc import Callable
from typing import Any, Literal

from .protocol import JsonObject

ApiMode = Literal['anthropic_messages', 'codex_responses', 'chat_completions']
# The values of the apiFormat container credential.
API_FORMATS: dict[str, ApiMode] = {'anthropic-messages': 'anthropic_messages', 'openai-responses': 'codex_responses',
                                   'openai-chat-completions': 'chat_completions'}
SCHEMA_NAME = 'output'
# Hermes sends no Responses call without instructions, and puts its own identity into empty ones.
NO_INSTRUCTIONS = 'You are a helpful assistant.'

# Keywords the strict structured-output grammar rejects; minLength and maxLength pass.
UNSUPPORTED = {
    'minimum', 'maximum', 'exclusiveMinimum', 'exclusiveMaximum', 'multipleOf',
    'maxItems', 'uniqueItems', 'contains', 'minContains', 'maxContains',
    'minProperties', 'maxProperties', 'patternProperties', 'propertyNames', 'dependentRequired',
    'dependentSchemas', 'unevaluatedProperties', 'unevaluatedItems', 'not', 'if', 'then', 'else',
}
STRING_FORMATS = {'date-time', 'time', 'date', 'duration', 'email', 'hostname', 'uri', 'ipv4', 'ipv6', 'uuid'}


def strict_schema(schema: Any) -> Any:
    """The schema in the strict grammar's subset: objects closed, `oneOf` as `anyOf`, unsupported keywords and
    formats dropped, `minItems` only 0 or 1. The reply is still checked against the original."""
    if isinstance(schema, list):
        return [strict_schema(s) for s in schema]
    if not isinstance(schema, dict):
        return schema
    out: JsonObject = {}
    for k, v in schema.items():
        if k in UNSUPPORTED or k == 'additionalProperties' or (k == 'format' and v not in STRING_FORMATS) or \
                (k == 'minItems' and isinstance(v, int) and v > 1):
            continue
        if k in ('properties', '$defs', 'definitions') and isinstance(v, dict):
            out[k] = {name: strict_schema(s) for name, s in v.items()}
        elif k in ('anyOf', 'allOf', 'prefixItems', 'items', 'oneOf'):
            out['anyOf' if k == 'oneOf' else k] = strict_schema(v)
        else:
            out[k] = v
    t = out.get('type')
    if t == 'object' or (isinstance(t, list) and 'object' in t) or (t is None and 'properties' in out):
        out['additionalProperties'] = False
    return out


def _anthropic_messages(request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
    out = dict(request)
    if system_prompt:
        out['system'] = [{'type': 'text', 'text': system_prompt, 'cache_control': {'type': 'ephemeral'}}]
    else:
        out.pop('system', None)
    if output_schema is not None:
        out['output_config'] = {**(out.get('output_config') or {}),
                                'format': {'type': 'json_schema', 'schema': strict_schema(output_schema)}}
    return out


def _openai_schema(output_schema: JsonObject) -> JsonObject:
    # Not strict: OpenAI's strict mode needs every property required; the reply is validated and retried instead.
    return {'name': SCHEMA_NAME, 'schema': output_schema, 'strict': False}


def _responses(request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
    # Nothing is stored: every call carries the whole history, its reasoning included.
    out = {**request, 'store': False, 'instructions': system_prompt or NO_INSTRUCTIONS}
    if output_schema is not None:
        out['text'] = {**(out.get('text') or {}), 'format': {'type': 'json_schema', **_openai_schema(output_schema)}}
    return out


def _chat_completions(request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
    out = dict(request)
    messages = list(out.get('messages') or [])
    # Hermes leads with a system message, or a developer one for the models that take that role.
    lead = messages.pop(0) if messages and messages[0].get('role') in ('system', 'developer') else {'role': 'system'}
    out['messages'] = ([{'role': lead['role'], 'content': system_prompt}] if system_prompt else []) + messages
    if output_schema is not None:
        out['response_format'] = {'type': 'json_schema', 'json_schema': _openai_schema(output_schema)}
    return out


PATCHES: dict[str, Callable[[JsonObject, str, JsonObject | None], JsonObject]] = {
    'anthropic_messages': _anthropic_messages, 'codex_responses': _responses, 'chat_completions': _chat_completions}
_REFUSED_MESSAGES = {'model': 'grokky-refused', 'max_tokens': 1, 'messages': [{'role': 'user', 'content': 'refused'}]}
REFUSED: dict[str, JsonObject] = {
    'anthropic_messages': _REFUSED_MESSAGES, 'chat_completions': _REFUSED_MESSAGES,
    'codex_responses': {'model': 'grokky-refused', 'instructions': 'refused', 'max_output_tokens': 16, 'store': False,
                        'input': [{'role': 'user', 'content': 'refused'}]}}


def patch(api_mode: str, request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
    """The outgoing request with the turn's system prompt, if any, and output format."""
    return PATCHES[api_mode](request, system_prompt, output_schema)


def refused(api_mode: str) -> JsonObject:
    """Sent instead of a request the session cannot vouch for: the API rejects it, and nothing of the conversation
    leaves the process. Middleware that raises lets Hermes send its own request, so an unknown format gets one too."""
    return copy.deepcopy(REFUSED.get(api_mode, REFUSED['anthropic_messages']))
