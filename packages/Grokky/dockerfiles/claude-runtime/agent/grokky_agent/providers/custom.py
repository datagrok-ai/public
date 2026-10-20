import re
from typing import ClassVar

from .. import models
from ..formats import API_FORMATS, ApiMode
from ..hermes import Connection
from .base import Credentials, Provider, ProviderError

# The model clients append their route to the path and cannot carry user info, a query or a fragment. A host with a
# trailing dot or outside ASCII is refused: Hermes drops the dot, and IDNA may map a fullwidth letter onto its ASCII
# one, either of which can turn another host into api.anthropic.com. Only https: over plain http, to a local host too,
# a tool of any session could listen on that port and take the key and the conversations.
BASE_URL = re.compile(r'https://[\w-]+(\.[\w-]+)*(:\d+)?(/[^?#\s]*)?', re.ASCII)
# What to check after a failed call, by its HTTP status; without one nothing answered as the API (a connection error,
# a login page instead of an event stream).
HINTS: dict[int | None, str] = {
    None: 'the baseUrl container credential',
    401: 'the apiKey container credential: a key the endpoint takes as a Bearer token',
    404: 'the baseUrl, apiFormat and model container credentials'}


class CustomProvider(Provider):
    """The customer's own gateway or self-hosted service at `baseUrl`, called in the `apiFormat` the admin names for
    every model, with the key as a Bearer token."""

    name = 'Custom'
    # Any catalog model by its own id, or a model the admin maps as `<id>=<catalog model>`; only the `model` and
    # `fastModel` credentials' models are offered.
    model_ids: ClassVar[dict[str, str]] = {m: m for m in models.CATALOG}

    def __init__(self, credentials: Credentials) -> None:
        super().__init__(credentials)
        self.models = {m: self.models[m] for m in (self.default, self.fast)}
        api_format = credentials.get('apiFormat') or ''
        if api_format not in API_FORMATS:
            raise ProviderError(f'The apiFormat container credential must be one of {", ".join(API_FORMATS)}')
        self._api_mode = API_FORMATS[api_format]
        url = credentials.get('baseUrl') or ''
        if not BASE_URL.fullmatch(url):
            raise ProviderError("The baseUrl container credential must be the endpoint's https URL, without user "
                                'info, a query or a fragment')
        self._base_url = self.endpoint(url)
        self.refuse_subscription_token(credentials.get('apiKey') or '', self._base_url)
        self._key = self.api_key(credentials, 'the key the endpoint takes as a Bearer token')

    def api_mode(self, model: str) -> ApiMode:
        return self._api_mode

    def connection(self, model: str) -> Connection:
        # On the Anthropic format Hermes sends a callable's key as Bearer, a string's as x-api-key; to
        # api.anthropic.com itself it sends a Console key as x-api-key either way.
        return {'provider': 'custom', 'api_mode': self._api_mode, 'base_url': self._base_url,
                'api_key': (lambda: self._key) if self._api_mode == 'anthropic_messages' else self._key}

    def error_hint(self, status: int | None) -> str:
        return f' (check {HINTS[status]})' if status in HINTS else ''
