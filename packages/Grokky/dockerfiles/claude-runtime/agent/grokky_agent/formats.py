"""The wire formats a model is called in, one class each, named by the api_mode Hermes uses for them: where a turn's
system prompt and output schema go on a request, and the request sent instead of one the session cannot vouch for."""
import copy
from abc import ABC, abstractmethod
from typing import Any, ClassVar, Literal

from .protocol import JsonObject

ApiMode = Literal['anthropic_messages', 'codex_responses', 'chat_completions', 'bedrock_converse']
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


def _openai_schema(output_schema: JsonObject) -> JsonObject:
    # Not strict: OpenAI's strict mode needs every property required; the reply is validated and retried instead.
    return {'name': SCHEMA_NAME, 'schema': output_schema, 'strict': False}


class ApiFormat(ABC):
    mode: ClassVar[ApiMode]
    # The value of the apiFormat container credential that picks this format for a custom endpoint.
    credential: ClassVar[str]
    refusal: ClassVar[JsonObject]

    @abstractmethod
    def patch(self, request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
        """The outgoing request with the turn's system prompt, if any, and output format."""

    def refused(self) -> JsonObject:
        """Sent instead of a request the session cannot vouch for: the API rejects it, and nothing of the
        conversation leaves the process."""
        return copy.deepcopy(self.refusal)


_REFUSED_MESSAGES: JsonObject = {'model': 'grokky-refused', 'max_tokens': 1,
                                 'messages': [{'role': 'user', 'content': 'refused'}]}


class AnthropicMessages(ApiFormat):
    mode = 'anthropic_messages'
    credential = 'anthropic-messages'
    refusal: ClassVar[JsonObject] = _REFUSED_MESSAGES

    def patch(self, request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
        out = dict(request)
        if system_prompt:
            out['system'] = [{'type': 'text', 'text': system_prompt, 'cache_control': {'type': 'ephemeral'}}]
        else:
            out.pop('system', None)
        if output_schema is not None:
            out['output_config'] = {**(out.get('output_config') or {}),
                                    'format': {'type': 'json_schema', 'schema': strict_schema(output_schema)}}
        return out


class OpenAIResponses(ApiFormat):
    mode = 'codex_responses'
    credential = 'openai-responses'
    refusal: ClassVar[JsonObject] = {'model': 'grokky-refused', 'instructions': 'refused', 'max_output_tokens': 16,
                                     'store': False, 'input': [{'role': 'user', 'content': 'refused'}]}

    def patch(self, request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
        # Nothing is stored: every call carries the whole history, its reasoning included.
        out = {**request, 'store': False, 'instructions': system_prompt or NO_INSTRUCTIONS}
        if output_schema is not None:
            out['text'] = {**(out.get('text') or {}),
                           'format': {'type': 'json_schema', **_openai_schema(output_schema)}}
        # Grokky shows no reasoning summary, and summaries may need a verified organization. Some models always reason
        # and reject `none`, which Hermes asks for when a turn wants no thinking; every reasoning model takes `low`.
        reasoning = {k: v for k, v in (out.pop('reasoning', None) or {}).items() if k != 'summary'}
        if reasoning.get('effort') == 'none':
            reasoning['effort'] = 'low'
        if reasoning:
            out['reasoning'] = reasoning
        return out


class ChatCompletions(ApiFormat):
    mode = 'chat_completions'
    credential = 'openai-chat-completions'
    refusal: ClassVar[JsonObject] = _REFUSED_MESSAGES

    def patch(self, request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
        out = dict(request)
        messages = list(out.get('messages') or [])
        # Hermes leads with a system message, or a developer one for the models that take that role.
        lead = messages.pop(0) if messages and messages[0].get('role') in ('system', 'developer') \
            else {'role': 'system'}
        out['messages'] = ([{'role': lead['role'], 'content': system_prompt}] if system_prompt else []) + messages
        if output_schema is not None:
            out['response_format'] = {'type': 'json_schema', 'json_schema': _openai_schema(output_schema)}
        return out


class BedrockConverse(ApiFormat):
    """Amazon Bedrock's Converse API, which parses the tool calls of models whose OpenAI-compatible hosts return them
    as text (Llama). Hermes calls it through boto3 with Converse's own request shape."""
    mode = 'bedrock_converse'
    credential = ''
    refusal: ClassVar[JsonObject] = {'modelId': 'grokky-refused', 'inferenceConfig': {'maxTokens': 1},
                                     'messages': [{'role': 'user', 'content': [{'text': 'refused'}]}]}

    def patch(self, request: JsonObject, system_prompt: str, output_schema: JsonObject | None) -> JsonObject:
        out = dict(request)
        if system_prompt:
            out['system'] = [{'text': system_prompt}]
        else:
            out.pop('system', None)
        return out


FORMATS: dict[str, ApiFormat] = {f.mode: f for f in (AnthropicMessages(), OpenAIResponses(), ChatCompletions(),
                                                    BedrockConverse())}
# The apiFormat container credential's values.
BY_CREDENTIAL: dict[str, ApiFormat] = {f.credential: f for f in FORMATS.values() if f.credential}


def refused(api_mode: str) -> JsonObject:
    """The refusal in the request's format; middleware that raises lets Hermes send its own request, so an unknown
    format gets one too."""
    return FORMATS.get(api_mode, FORMATS['anthropic_messages']).refused()
