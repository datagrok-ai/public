from collections.abc import Callable

from .anthropic import AnthropicProvider
from .base import Provider, ProviderError

PROVIDERS: dict[str, Callable[[dict[str, str]], Provider]] = {'Anthropic': AnthropicProvider}


def from_credentials(credentials: dict[str, str]) -> Provider:
    name = credentials.get('provider') or 'Anthropic'
    if name not in PROVIDERS:
        raise ProviderError(f'Unsupported provider "{name}" in the container credentials')
    return PROVIDERS[name](credentials)


__all__ = ['Provider', 'ProviderError', 'from_credentials']
