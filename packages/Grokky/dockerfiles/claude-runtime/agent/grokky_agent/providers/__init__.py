"""The providers the container credentials can name: each module of this package defines one Provider subclass, which
registers it under its `name`."""
import importlib
import pkgutil

from .base import Credentials, Provider, ProviderError

for _module in pkgutil.iter_modules(__path__):
    importlib.import_module(f'{__name__}.{_module.name}')

PROVIDERS: dict[str, type[Provider]] = {p.name: p for p in Provider.__subclasses__()}


def from_credentials(credentials: Credentials) -> Provider:
    name = credentials.get('provider') or 'Anthropic'
    if name not in PROVIDERS:
        raise ProviderError(f'Unsupported provider "{name}" in the container credentials')
    return PROVIDERS[name](credentials)


__all__ = ['PROVIDERS', 'Provider', 'ProviderError', 'from_credentials']
