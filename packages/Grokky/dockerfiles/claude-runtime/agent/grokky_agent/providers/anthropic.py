from ..formats import ApiMode
from .base import Credentials, Provider, Route

API = 'https://api.anthropic.com'


class AnthropicProvider(Provider):
    """Claude on the Anthropic Messages API with a Console API key."""

    name = 'Anthropic'
    key_kind = 'an Anthropic Console API key'
    model_example = 'claude-sonnet-5-5'
    billing = 'anthropic'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        url = self.endpoint(API)
        self.refuse_subscription_token(credentials.get('apiKey') or '', url)
        return {'anthropic_messages': Route(url)}

    def error_hint(self, status: int | None) -> str:
        return ' (check the apiKey container credential: an Anthropic Console API key)' if status == 401 else ''
