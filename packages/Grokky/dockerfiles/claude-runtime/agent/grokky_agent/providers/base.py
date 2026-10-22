import re
from abc import ABC, abstractmethod
from typing import ClassVar, NamedTuple
from urllib.parse import urlsplit, urlunsplit

from .. import hermes
from ..formats import FORMATS, ApiFormat, ApiMode
from ..hermes import Connection, Usage
from ..protocol import JsonObject

Credentials = dict[str, str]

# Claude speaks Anthropic's Messages API everywhere; GPT calls tools reliably only through OpenAI's Responses API.
CLAUDE = re.compile(r'claude', re.IGNORECASE)
LOCAL_HOST = r'(127\.0\.0\.1|localhost|\[::1\]|host\.docker\.internal)'
# Hermes sends any other sk-ant-* key, a JWT or a cc- token to api.anthropic.com as a Claude subscription login under
# Claude Code's identity (agent/anthropic_credentials.py _is_oauth_token), which Anthropic allows only in its own apps.
SUBSCRIPTION_TOKEN = re.compile(r'^(sk-ant-(?!api|usr|svc|admin)|eyJ|cc-)')


class ProviderError(ValueError):
    """The container credentials cannot be used with this provider."""


class Route(NamedTuple):
    """Where a provider serves a wire format, and whether the key goes out as a Bearer token. On the Anthropic
    Messages format Hermes sends a string key as x-api-key and a callable's as Bearer (hermes.send_keys_by_type);
    the OpenAI formats always send Bearer."""
    url: str
    bearer: bool = False


class Provider(ABC):
    """A model API the agent runs on: the container credentials it takes, the wire formats it serves and where, and
    the one model every turn runs on, named by the `model` credential as the provider names it."""

    name: ClassVar[str]
    # What the apiKey credential must be, for the message when it is missing.
    key_kind: ClassVar[str]
    # A model as this provider names it, for the message when the model credential is missing.
    model_example: ClassVar[str]
    # The price list Hermes keeps for this provider (agent/usage_pricing.py). Without one, or for a model the list
    # lacks, a turn's cost is unknown and only its call limit bounds it.
    billing: ClassVar[str | None] = None

    def __init__(self, credentials: Credentials) -> None:
        self.test_upstream = credentials.get('TEST_UPSTREAM')
        self.model = self._model(credentials)
        routes = self.routes(credentials)
        mode = next(iter(routes)) if len(routes) == 1 else \
            'anthropic_messages' if CLAUDE.search(self.model) else 'codex_responses'
        self.format: ApiFormat = FORMATS[mode]
        self.route = routes[mode]
        key = credentials.get('apiKey')
        if not key:
            raise ProviderError(f'The apiKey container credential is not set: set {self.key_kind}')
        if re.search(r'[^\x21-\x7e]', key):
            raise ProviderError('The apiKey container credential contains whitespace or control characters')
        self.check_key(key)
        self.key = key

    @abstractmethod
    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        """The wire formats this provider serves, each with where; with two, a Claude model is called in Anthropic's
        Messages format and any other in OpenAI's Responses one."""

    def check_key(self, key: str) -> None:
        """Refuses a key of the wrong kind before a call is made."""
        return None

    def connection(self) -> Connection:
        """How Hermes reaches the model: always an explicit Hermes provider id, never `auto` or empty, with which
        Hermes would pick a provider and a credential itself."""
        key = self.key
        bearer = self.route.bearer and self.format.mode == 'anthropic_messages'
        return {'provider': 'custom', 'api_mode': self.format.mode, 'base_url': self.route.url,
                'api_key': (lambda: key) if bearer else key}

    def patch_request(self, api_mode: str, request: JsonObject, system_prompt: str,
                      output_schema: JsonObject | None) -> JsonObject:
        """The outgoing request with the turn's system prompt, if any, and output format."""
        return FORMATS[api_mode].patch(request, system_prompt, output_schema)

    def cost(self, usage: Usage) -> float | None:
        """USD for one call, or None when Hermes has no price for the model."""
        return hermes.usage_cost(self.model, usage, self.billing, self.route.url) if self.billing else None

    def error_hint(self, status: int | None) -> str:
        """Added to a failed call's message, e.g. what to check for a rejected key."""
        return ''

    def endpoint(self, url: str) -> str:
        """The URL, or with TEST_UPSTREAM set, the URL with its scheme, host and port replaced: the path stays, so a
        test stub sees each route's real path. The override may name only a local host, or this URL's own host over
        https for a stub answering under that name."""
        if not self.test_upstream:
            return url
        real = urlsplit(url)
        if not re.fullmatch(rf'(https?://{LOCAL_HOST}|https://{re.escape(real.hostname or "")})(:\d+)?/?',
                            self.test_upstream):
            raise ProviderError("TEST_UPSTREAM may only name a local test host, or the provider's own host over https")
        upstream = urlsplit(self.test_upstream)
        return urlunsplit(real._replace(scheme=upstream.scheme, netloc=upstream.netloc))

    @staticmethod
    def refuse_subscription_token(key: str, url: str) -> None:
        """Refuses a Claude subscription token bound for api.anthropic.com; other hosts take JWTs (Entra ID, OAuth)."""
        if urlsplit(url).hostname == 'api.anthropic.com' and SUBSCRIPTION_TOKEN.match(key.strip()):
            raise ProviderError('The apiKey container credential is a Claude subscription token: '
                                'set an Anthropic Console API key')

    def _model(self, credentials: Credentials) -> str:
        model = (credentials.get('model') or '').strip()
        if not model:
            raise ProviderError(f'The model container credential is not set: name the model as {self.name} names '
                                f'it, such as {self.model_example}')
        if re.search(r'\s', model):
            raise ProviderError('The model container credential contains whitespace')
        return model
