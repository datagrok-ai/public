from typing import ClassVar

from ..hermes import Connection
from ..protocol import JsonObject
from .base import Credentials, Provider, ProviderError

API = 'https://openrouter.ai/api'
# Replaces Hermes' own app attribution (https://openrouter.ai/docs/app-attribution); an empty category list drops the
# categories Hermes files its calls under.
ATTRIBUTION = {'HTTP-Referer': 'https://datagrok.ai', 'X-Title': 'Datagrok', 'X-OpenRouter-Categories': ''}


class OpenRouterProvider(Provider):
    """Claude on OpenRouter's Anthropic Messages API and GPT on its Responses API, with an OpenRouter API key."""

    name = 'OpenRouter'
    # https://openrouter.ai/api/v1/models
    model_ids: ClassVar[dict[str, str]] = {
        'anthropic/claude-haiku-5.5': 'claude-haiku-5-5', 'anthropic/claude-sonnet-5.5': 'claude-sonnet-5-5',
        'anthropic/claude-opus-5.5': 'claude-opus-5-5', 'anthropic/claude-fable-5.1': 'claude-fable-5-1',
        **{f'openai/{m}': m for m in ('gpt-6-luna', 'gpt-6-sol', 'gpt-6.1-sol', 'gpt-6-astra', 'gpt-5.6-luna',
                                      'gpt-5.6-terra', 'gpt-5.6-sol')}}
    default_model, fast_model = 'anthropic/claude-sonnet-5.5', 'anthropic/claude-haiku-5.5'

    def __init__(self, credentials: Credentials) -> None:
        super().__init__(credentials)
        self._messages_url = self.endpoint(API)
        self._responses_url = self.endpoint(f'{API}/v1')
        self._key = self.api_key(credentials, 'an OpenRouter API key')
        if not self._key.startswith('sk-or-'):
            raise ProviderError('The apiKey container credential is not an OpenRouter API key (sk-or-…)')
        # Hermes names Claude on the Messages API by Anthropic's own id (claude-sonnet-5-5 for
        # anthropic/claude-sonnet-5.5); the request goes out with the OpenRouter id it runs.
        self._ids = {m.removeprefix('anthropic/').replace('.', '-'): m
                     for m in self.models if self.api_mode(m) == 'anthropic_messages'}

    def connection(self, model: str) -> Connection:
        if self.api_mode(model) == 'anthropic_messages':
            # A callable key goes out as Bearer, OpenRouter's documented scheme.
            return {'provider': 'custom', 'api_mode': 'anthropic_messages', 'base_url': self._messages_url,
                    'api_key': lambda: self._key}
        return {'provider': 'custom', 'api_mode': 'codex_responses', 'base_url': self._responses_url,
                'api_key': self._key}

    def patch_request(self, api_mode: str, request: JsonObject, system_prompt: str,
                      output_schema: JsonObject | None) -> JsonObject:
        out = super().patch_request(api_mode, request, system_prompt, output_schema)
        out['extra_headers'] = {**(out.get('extra_headers') or {}), **ATTRIBUTION}
        if api_mode == 'anthropic_messages':
            out['model'] = self._ids.get(out['model'], out['model'])
            # Some hosts of every Claude model lack structured output; GPT's lack it only where they lack
            # response_format, which OpenRouter already routes around.
            if output_schema is not None:
                out['extra_body'] = {**(out.get('extra_body') or {}), 'provider': {'require_parameters': True}}
        return out

    def error_hint(self, status: int | None) -> str:
        return ' (check the apiKey container credential: an OpenRouter API key)' if status == 401 else ''
