import re

from ..formats import ApiMode
from .base import CLAUDE, Credentials, Provider, ProviderError, Route

# A Microsoft Foundry resource serves Claude and GPT, an Azure OpenAI resource GPT only. Claude speaks the Messages API
# under /anthropic, GPT the v1 API under /openai/v1, which takes no api-version:
# https://learn.microsoft.com/en-us/azure/foundry/foundry-models/how-to/use-foundry-models-claude
# https://learn.microsoft.com/en-us/azure/foundry/openai/api-version-lifecycle
RESOURCE = re.compile(r'https://[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\.(services\.ai|openai)\.azure\.com')
# Azure answers an unknown deployment with 404.
HINTS: dict[int | None, str] = {
    401: ' (check the apiKey container credential: a key of the Azure resource in baseUrl)',
    404: ' (check that baseUrl is the Azure resource and that model names one of its deployments)'}


class AzureProvider(Provider):
    """Claude on Microsoft Foundry and GPT on Azure OpenAI, called by deployment name with the resource's key."""

    name = 'Azure'
    key_kind = 'a key of the Azure resource'
    model_example = 'claude-sonnet-5-5 (a deployment name)'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        resource = RESOURCE.fullmatch((credentials.get('baseUrl') or '').strip().rstrip('/').lower())
        if not resource:
            raise ProviderError('The baseUrl container credential must be the Azure resource: '
                                'https://<resource>.services.ai.azure.com, or https://<resource>.openai.azure.com '
                                'for GPT only')
        if CLAUDE.search(self.model):
            if resource[1] == 'openai':
                raise ProviderError('The model container credential names a Claude deployment, which needs a '
                                    'Microsoft Foundry resource, https://<resource>.services.ai.azure.com: an Azure '
                                    'OpenAI resource serves GPT only')
            # Hermes picks a Claude request's thinking and sampling settings by the model name, and turns the dots
            # of a claude-* name into hyphens.
            if '.' in self.model:
                raise ProviderError('The model container credential names a Claude deployment with a dot in its '
                                    'name: name the deployment after its model id, as Azure does by default')
        routes: dict[ApiMode, Route] = {'codex_responses': Route(self.endpoint(resource[0] + '/openai/v1'))}
        if resource[1] == 'services.ai':
            routes['anthropic_messages'] = Route(self.endpoint(resource[0] + '/anthropic'))
        return routes

    def error_hint(self, status: int | None) -> str:
        return HINTS.get(status, '')
