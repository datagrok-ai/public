import re

from .. import models
from ..formats import ApiMode
from ..hermes import Connection
from .base import Credentials, Provider, ProviderError

# A Microsoft Foundry resource serves Claude and GPT, an Azure OpenAI resource GPT only. Claude speaks the Messages API
# under /anthropic, GPT the v1 API under /openai/v1, which takes no api-version:
# https://learn.microsoft.com/en-us/azure/foundry/foundry-models/how-to/use-foundry-models-claude
# https://learn.microsoft.com/en-us/azure/foundry/openai/api-version-lifecycle
RESOURCE = re.compile(r'https://[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\.(services\.ai|openai)\.azure\.com')
ROUTES: dict[ApiMode, str] = {'anthropic_messages': '/anthropic', 'codex_responses': '/openai/v1'}
# Azure answers an unknown deployment with 404.
HINTS: dict[int | None, str] = {
    401: ' (check the apiKey container credential: a key of the Azure resource in baseUrl)',
    404: ' (check that baseUrl is the Azure resource and that model and fastModel name its deployments)'}


class AzureProvider(Provider):
    """Claude on Microsoft Foundry and GPT on Azure OpenAI, called by deployment name with the resource's key."""

    name = 'Azure'

    def __init__(self, credentials: Credentials) -> None:
        # A deployment is named after its model by default: a bare catalog id names a deployment of that model.
        super().__init__({**credentials, **{c: f'{v}={v}' for c in ('model', 'fastModel')
                                            if (v := (credentials.get(c) or '').strip()) in models.CATALOG}})
        resource = RESOURCE.fullmatch((credentials.get('baseUrl') or '').strip().rstrip('/').lower())
        if not resource:
            raise ProviderError('The baseUrl container credential must be the Azure resource: '
                                'https://<resource>.services.ai.azure.com, or https://<resource>.openai.azure.com '
                                'for GPT only')
        claude = sorted(m for m, name in self.models.items() if models.CATALOG[name]['family'] == 'claude')
        if claude and resource[1] == 'openai':
            raise ProviderError(f'Claude ({", ".join(claude)}) needs a Microsoft Foundry resource, '
                                'https://<resource>.services.ai.azure.com: an Azure OpenAI resource serves GPT only')
        # Hermes picks a Claude request's thinking and sampling settings by the model name, and turns the dots of a
        # claude-* name into hyphens: under another name Opus 5.5 and Fable 5.1 would get thinking they refuse.
        for d in claude:
            if self.models[d] not in d or '.' in d:
                raise ProviderError(f'The Claude deployment "{d}" needs a name that contains its model id and no dot, '
                                    f'as Azure names it by default ({self.models[d]})')
        self._base_urls = {mode: self.endpoint(resource[0] + path) for mode, path in ROUTES.items()}
        self._key = self.api_key(credentials, 'a key of the Azure resource')

    def connection(self, model: str) -> Connection:
        mode = self.api_mode(model)
        return {'provider': 'custom', 'api_mode': mode, 'base_url': self._base_urls[mode], 'api_key': self._key}

    def error_hint(self, status: int | None) -> str:
        return HINTS.get(status, '')
