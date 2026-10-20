import re
from abc import ABC, abstractmethod
from typing import ClassVar
from urllib.parse import urlsplit, urlunsplit

from .. import formats, models
from ..formats import ApiMode
from ..hermes import Connection
from ..protocol import JsonObject

Credentials = dict[str, str]

# Claude speaks Anthropic's Messages API everywhere; GPT calls tools reliably only through OpenAI's Responses API.
FAMILY_API_MODES: dict[str, ApiMode] = {'claude': 'anthropic_messages', 'gpt': 'codex_responses'}
LOCAL_HOST = r'(127\.0\.0\.1|localhost|\[::1\]|host\.docker\.internal)'
# Hermes sends any other sk-ant-* key, a JWT or a cc- token to api.anthropic.com as a Claude subscription login under
# Claude Code's identity (agent/anthropic_credentials.py _is_oauth_token), which Anthropic allows only in its own apps.
SUBSCRIPTION_TOKEN = re.compile(r'^(sk-ant-(?!api|usr|svc|admin)|eyJ|cc-)')


class ProviderError(ValueError):
    """The container credentials cannot be used with this provider."""


class Provider(ABC):
    """A model API the agent runs on: it checks the container credentials, names its models and the catalog model
    each one runs (models.yaml), and tells Hermes how to reach each model."""

    name: ClassVar[str]
    # The provider's model ids, each with the catalog model it runs, and its default and fast model. The `model` and
    # `fastModel` credentials change the two, and add a model the admin names as `<id>=<catalog model>`; `model`
    # alone is the fast model too, as the provider's own may run in another region or at another vendor.
    model_ids: ClassVar[dict[str, str]] = {}
    default_model: ClassVar[str | None] = None
    fast_model: ClassVar[str | None] = None

    def __init__(self, credentials: Credentials) -> None:
        self.test_upstream = credentials.get('TEST_UPSTREAM')
        self.models = dict(self.model_ids)
        self.default = self._model(credentials, 'model', self.default_model)
        self.fast = self._model(credentials, 'fastModel',
                                self.fast_model if self.fast_model and not credentials.get('model') else self.default)

    @abstractmethod
    def connection(self, model: str) -> Connection:
        """How Hermes reaches the model: always an explicit Hermes provider id, never `auto` or empty, with which
        Hermes would pick a provider and a credential itself."""

    def api_mode(self, model: str) -> ApiMode:
        """The request format the model is called in."""
        return FAMILY_API_MODES[models.CATALOG[self.models[model]]['family']]

    def images(self, model: str) -> bool:
        """Whether the model takes image input on this provider."""
        return models.CATALOG[self.models[model]]['images']

    def patch_request(self, api_mode: str, request: JsonObject, system_prompt: str,
                      output_schema: JsonObject | None) -> JsonObject:
        """The outgoing request with the turn's system prompt, if any, and output format."""
        out = formats.patch(api_mode, request, system_prompt, output_schema)
        if api_mode != 'codex_responses':
            return out
        # Grokky shows no reasoning summary, and summaries may need a verified organization. Hermes asks a model it
        # does not know by name (a provider's prefixed id, a deployment) for `effort: none`, which one that always
        # reasons rejects; without an effort it reasons at its default.
        reasoning = {k: v for k, v in (out.pop('reasoning', None) or {}).items() if k != 'summary'}
        catalog = models.CATALOG.get(self.models.get(out.get('model', ''), ''))
        if reasoning.get('effort') == 'none' and catalog and catalog.get('always_reasons'):
            reasoning = {}
        if reasoning:
            out['reasoning'] = reasoning
        return out

    def refused_request(self, api_mode: str) -> JsonObject:
        return formats.refused(api_mode)

    def error_hint(self, status: int | None) -> str:
        """Added to a failed call's message, e.g. what to check for a rejected key."""
        return ''

    def hermes_models(self) -> JsonObject:
        """Hermes' metadata of each model (config.yaml `model_overrides`), by Hermes provider id."""
        overrides: dict[str, JsonObject] = {}
        for model, name in self.models.items():
            overrides.setdefault(self.connection(model)['provider'], {})[model] = {
                'context_window': models.CATALOG[name]['context'], 'supports_vision': self.images(model)}
        return overrides

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
    def api_key(credentials: Credentials, what: str) -> str:
        key = credentials.get('apiKey')
        if not key:
            raise ProviderError(f'The apiKey container credential is not set: set {what}')
        if re.search(r'[^\x21-\x7e]', key):
            raise ProviderError('The apiKey container credential contains whitespace or control characters')
        return key

    @staticmethod
    def refuse_subscription_token(key: str, url: str) -> None:
        """Refuses a Claude subscription token bound for api.anthropic.com; other hosts take JWTs (Entra ID, OAuth)."""
        if urlsplit(url).hostname == 'api.anthropic.com' and SUBSCRIPTION_TOKEN.match(key.strip()):
            raise ProviderError('The apiKey container credential is a Claude subscription token: '
                                'set an Anthropic Console API key')

    def _model(self, credentials: Credentials, credential: str, fallback: str | None) -> str:
        """The model a credential names, as one of this provider's ids or `id=catalog model`, else `fallback`."""
        value = credentials.get(credential)
        if not value:
            if fallback is None:
                raise ProviderError(f'The {credential} container credential is not set: '
                                    'name a model as <id>=<catalog model>')
            return fallback
        model, _, name = (part.strip() for part in value.partition('='))
        if name:
            if name not in models.CATALOG:
                raise ProviderError(f'The {credential} container credential must map an id to a catalog model, one of '
                                    f'{", ".join(models.CATALOG)}')
            self.models[model] = name
        elif model not in self.models:
            raise ProviderError(f'The {credential} container credential must be one of '
                                f'{", ".join([*self.models, "<id>=<catalog model>"])}')
        return model
