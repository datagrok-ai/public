from ..formats import ApiMode
from ..protocol import JsonObject
from .base import CLAUDE, Credentials, Provider, ProviderError, Route

API = 'https://openrouter.ai/api'
# Replaces Hermes' own app attribution (https://openrouter.ai/docs/app-attribution); an empty category list drops the
# categories Hermes files its calls under.
ATTRIBUTION = {'HTTP-Referer': 'https://datagrok.ai', 'X-Title': 'Datagrok', 'X-OpenRouter-Categories': ''}


class OpenRouterProvider(Provider):
    """Claude on OpenRouter's Anthropic Messages API, OpenAI's models on its Responses API and every other model on
    its Chat Completions API, with an OpenRouter API key, Bearer as OpenRouter documents it. Gemini, Qwen and Llama
    return no tool calls through the Responses API."""

    name = 'OpenRouter'
    key_kind = 'an OpenRouter API key'
    model_example = 'anthropic/claude-sonnet-5.5'
    billing = 'openrouter'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        return {'anthropic_messages': Route(self.endpoint(API), bearer=True),
                'codex_responses': Route(self.endpoint(f'{API}/v1')),
                'chat_completions': Route(self.endpoint(f'{API}/v1'))}

    def mode(self) -> ApiMode:
        if CLAUDE.search(self.model):
            return 'anthropic_messages'
        return 'codex_responses' if self.model.startswith('openai/') else 'chat_completions'

    def check_key(self, key: str) -> None:
        if not key.startswith('sk-or-'):
            raise ProviderError('The apiKey container credential is not an OpenRouter API key (sk-or-…)')

    def patch_request(self, api_mode: str, request: JsonObject, system_prompt: str,
                      output_schema: JsonObject | None) -> JsonObject:
        out = super().patch_request(api_mode, request, system_prompt, output_schema)
        out['extra_headers'] = {**(out.get('extra_headers') or {}), **ATTRIBUTION}
        # Some hosts of every Claude model lack structured output; GPT's lack it only where they lack
        # response_format, which OpenRouter already routes around.
        if api_mode == 'anthropic_messages' and output_schema is not None:
            out['extra_body'] = {**(out.get('extra_body') or {}), 'provider': {'require_parameters': True}}
        return out

    def error_hint(self, status: int | None) -> str:
        return ' (check the apiKey container credential: an OpenRouter API key)' if status == 401 else ''
