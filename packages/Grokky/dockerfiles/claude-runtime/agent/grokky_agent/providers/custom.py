import re

from ..formats import BY_CREDENTIAL, ApiMode
from .base import Credentials, Provider, ProviderError, Route

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
    """The customer's own gateway or self-hosted service at `baseUrl`, called in the `apiFormat` the admin names, with
    the key as a Bearer token."""

    name = 'Custom'
    key_kind = 'the key the endpoint takes as a Bearer token'
    model_example = 'claude-sonnet-5-5'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        api_format = BY_CREDENTIAL.get(credentials.get('apiFormat') or '')
        if api_format is None:
            raise ProviderError(f'The apiFormat container credential must be one of {", ".join(BY_CREDENTIAL)}')
        url = credentials.get('baseUrl') or ''
        if not BASE_URL.fullmatch(url):
            raise ProviderError("The baseUrl container credential must be the endpoint's https URL, without user "
                                'info, a query or a fragment')
        route = Route(self.endpoint(url), bearer=True)
        self.refuse_subscription_token(credentials.get('apiKey') or '', route.url)
        return {api_format.mode: route}

    def error_hint(self, status: int | None) -> str:
        return f' (check {HINTS[status]})' if status in HINTS else ''
