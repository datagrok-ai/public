from typing import ClassVar

from ..hermes import Connection
from .base import Credentials, Provider

API = 'https://api.anthropic.com'


class AnthropicProvider(Provider):
    """Claude on the Anthropic Messages API with a Console API key."""

    name = 'Anthropic'
    # https://platform.claude.com/docs/en/about-claude/models/overview
    model_ids: ClassVar[dict[str, str]] = {
        m: m for m in ('claude-haiku-5-5', 'claude-sonnet-5-5', 'claude-opus-5-5', 'claude-fable-5-1')}
    default_model, fast_model = 'claude-sonnet-5-5', 'claude-haiku-5-5'

    def __init__(self, credentials: Credentials) -> None:
        super().__init__(credentials)
        self._base_url = self.endpoint(API)
        self.refuse_subscription_token(credentials.get('apiKey') or '', self._base_url)
        self._key = self.api_key(credentials, 'an Anthropic Console API key')

    def connection(self, model: str) -> Connection:
        return {'provider': 'custom', 'api_mode': self.api_mode(model), 'base_url': self._base_url,
                'api_key': self._key}

    def error_hint(self, status: int | None) -> str:
        return ' (check the apiKey container credential: an Anthropic Console API key)' if status == 401 else ''
