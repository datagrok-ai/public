import re
from typing import Any

from ..hermes import Connection
from ..protocol import JsonObject
from .base import Provider, ProviderError

API = 'https://api.anthropic.com'
# Hermes sends any other sk-ant-* key, a JWT or a cc- token as a Claude subscription login under Claude Code's
# identity (agent/anthropic_credentials.py _is_oauth_token), which Anthropic allows only in its own applications.
SUBSCRIPTION_TOKEN = re.compile(r'^(sk-ant-(?!api|usr|svc|admin)|eyJ|cc-)')
# The upstream override exists for test stubs and loggers; plain http only to a local host.
TEST_UPSTREAM = re.compile(
    r'^(https://api\.anthropic\.com|https?://(127\.0\.0\.1|localhost|\[::1\]|host\.docker\.internal))(:\d+)?(/|$)')

# Keywords the strict structured-output grammar rejects; minLength and maxLength pass.
UNSUPPORTED = {
    'minimum', 'maximum', 'exclusiveMinimum', 'exclusiveMaximum', 'multipleOf',
    'maxItems', 'uniqueItems', 'contains', 'minContains', 'maxContains',
    'minProperties', 'maxProperties', 'patternProperties', 'propertyNames', 'dependentRequired',
    'dependentSchemas', 'unevaluatedProperties', 'unevaluatedItems', 'not', 'if', 'then', 'else',
}
FORMATS = {'date-time', 'time', 'date', 'duration', 'email', 'hostname', 'uri', 'ipv4', 'ipv6', 'uuid'}


def strict_schema(schema: Any) -> Any:
    """The schema in the strict grammar's subset: objects closed, `oneOf` as `anyOf`, unsupported keywords and
    formats dropped, `minItems` only 0 or 1. The reply is still checked against the original."""
    if isinstance(schema, list):
        return [strict_schema(s) for s in schema]
    if not isinstance(schema, dict):
        return schema
    out: JsonObject = {}
    for k, v in schema.items():
        if k in UNSUPPORTED or k == 'additionalProperties' or (k == 'format' and v not in FORMATS) or \
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


class AnthropicProvider(Provider):
    """The Anthropic Messages API with a Console API key."""

    def __init__(self, credentials: dict[str, str]) -> None:
        key = credentials.get('apiKey')
        if not key:
            raise ProviderError('The apiKey container credential is not set: set an Anthropic Console API key')
        if SUBSCRIPTION_TOKEN.match(key.strip()):
            raise ProviderError('The apiKey container credential is a Claude subscription token: '
                                'set an Anthropic Console API key')
        if re.search(r'[^\x21-\x7e]', key):
            raise ProviderError('The apiKey container credential contains whitespace or control characters')
        upstream = credentials.get('ANTHROPIC_UPSTREAM') or API
        if not TEST_UPSTREAM.match(upstream):
            raise ProviderError('ANTHROPIC_UPSTREAM may only name api.anthropic.com or a local test host')
        self._key = key
        self._base_url = upstream

    def connection(self) -> Connection:
        return {'provider': 'custom', 'api_mode': 'anthropic_messages', 'base_url': self._base_url,
                'api_key': self._key}

    def patch_request(self, request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
        out = dict(request)
        if system_prompt:
            out['system'] = [{'type': 'text', 'text': system_prompt, 'cache_control': {'type': 'ephemeral'}}]
        else:
            out.pop('system', None)
        if output_schema is not None:
            out['output_config'] = {**(out.get('output_config') or {}),
                                    'format': {'type': 'json_schema', 'schema': strict_schema(output_schema)}}
        return out

    def error_hint(self, status: int | None) -> str:
        return ' (check the apiKey container credential: an Anthropic Console API key)' if status == 401 else ''
