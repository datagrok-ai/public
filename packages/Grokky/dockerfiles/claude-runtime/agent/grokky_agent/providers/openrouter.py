from ..formats import ApiMode
from ..protocol import JsonObject
from .base import Credentials, Provider, ProviderError, Route

API = 'https://openrouter.ai/api'
# Replaces Hermes' own app attribution (https://openrouter.ai/docs/app-attribution); an empty category list drops the
# categories Hermes files its calls under.
ATTRIBUTION = {'HTTP-Referer': 'https://datagrok.ai', 'X-Title': 'Datagrok', 'X-OpenRouter-Categories': ''}


class OpenRouterProvider(Provider):
    """Claude on OpenRouter's Anthropic Messages API and any other model on its Responses API, with an OpenRouter
    API key, Bearer on both as OpenRouter documents it."""

    name = 'OpenRouter'
    key_kind = 'an OpenRouter API key'
    model_example = 'anthropic/claude-sonnet-5.5'
    billing = 'openrouter'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        return {'anthropic_messages': Route(self.endpoint(API), bearer=True),
                'codex_responses': Route(self.endpoint(f'{API}/v1'))}

    def check_key(self, key: str) -> None:
        if not key.startswith('sk-or-'):
            raise ProviderError('The apiKey container credential is not an OpenRouter API key (sk-or-…)')

    def patch_request(self, api_mode: str, request: JsonObject, system_prompt: str,
                      output_schema: JsonObject | None) -> JsonObject:
        out = super().patch_request(api_mode, request, system_prompt, output_schema)
        out['extra_headers'] = {**(out.get('extra_headers') or {}), **ATTRIBUTION}
        if api_mode == 'anthropic_messages':
            # Hermes names Claude by Anthropic's own id on the Messages API (claude-sonnet-5-5 for
            # anthropic/claude-sonnet-5.5); OpenRouter routes by its own.
            out['model'] = self.model
            # Some hosts of every Claude model lack structured output; GPT's lack it only where they lack
            # response_format, which OpenRouter already routes around.
            if output_schema is not None:
                out['extra_body'] = {**(out.get('extra_body') or {}), 'provider': {'require_parameters': True}}
        return out

    def error_hint(self, status: int | None) -> str:
        return ' (check the apiKey container credential: an OpenRouter API key)' if status == 401 else ''
