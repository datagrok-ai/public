import re
from typing import ClassVar

from ..hermes import Connection
from .base import Credentials, Provider, ProviderError

API = 'https://api.openai.com/v1'
# A project with data residency is called at its region's host:
# https://developers.openai.com/api/docs/guides/your-data
REGIONS = ('us', 'eu', 'ae', 'au', 'ca', 'jp', 'in', 'sg', 'kr', 'gb')
BASE_URL = re.compile(rf'https://((?:{"|".join(REGIONS)})\.)?api\.openai\.com(?:/v1)?/?')


class OpenAIProvider(Provider):
    """GPT on the OpenAI Responses API with a project's API key."""

    name = 'OpenAI'
    # https://developers.openai.com/api/docs/models/all
    model_ids: ClassVar[dict[str, str]] = {m: m for m in (
        'gpt-6-luna', 'gpt-6-sol', 'gpt-6.1-sol', 'gpt-6-astra', 'gpt-5.6-luna', 'gpt-5.6-terra', 'gpt-5.6-sol')}
    # OpenAI's balance of intelligence and cost, near Astra at a fifth of its price (Astra may not be open to every
    # account), and its cheapest model.
    default_model, fast_model = 'gpt-6.1-sol', 'gpt-6-luna'

    def __init__(self, credentials: Credentials) -> None:
        super().__init__(credentials)
        self._base_url = self.endpoint(self._api_url(credentials.get('baseUrl') or ''))
        key = (credentials.get('apiKey') or '').strip()
        # A ChatGPT or Codex login and a workload identity federation token are JWTs.
        if key.startswith('eyJ'):
            raise ProviderError('The apiKey container credential is a token, such as a ChatGPT or Codex login: '
                                'set an API key of an OpenAI project')
        if key.startswith('sk-admin-'):
            raise ProviderError('The apiKey container credential is an OpenAI Admin API key, which works only on the '
                                'administration endpoints: set an API key of a project')
        self._key = self.api_key(credentials, 'an API key of an OpenAI project')

    @staticmethod
    def _api_url(base_url: str) -> str:
        if not base_url:
            return API
        m = BASE_URL.fullmatch(base_url.strip())
        if not m:
            raise ProviderError('The baseUrl container credential is not an OpenAI API address: set '
                                f'https://<region>.api.openai.com, region one of {", ".join(REGIONS)}, for a project '
                                'with data residency, or leave it empty')
        return f'https://{m[1] or ""}api.openai.com/v1'

    def connection(self, model: str) -> Connection:
        return {'provider': 'custom', 'api_mode': self.api_mode(model), 'base_url': self._base_url,
                'api_key': self._key}

    def error_hint(self, status: int | None) -> str:
        return (' (check the apiKey container credential, an API key of an OpenAI project, and for a project with '
                'data residency the baseUrl one)') if status == 401 else ''
