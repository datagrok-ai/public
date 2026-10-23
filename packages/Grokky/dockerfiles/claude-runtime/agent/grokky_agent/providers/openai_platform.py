import re

from ..formats import ApiMode
from .base import Credentials, Provider, ProviderError, Route

API = 'https://api.openai.com/v1'
# A project with data residency is called at its region's host:
# https://developers.openai.com/api/docs/guides/your-data
REGIONS = ('us', 'eu', 'ae', 'au', 'ca', 'jp', 'in', 'sg', 'kr', 'gb')
BASE_URL = re.compile(rf'https://((?:{"|".join(REGIONS)})\.)?api\.openai\.com(?:/v1)?/?')


class OpenAIProvider(Provider):
    """GPT on the OpenAI Responses API with a project's API key."""

    name = 'OpenAI'
    key_kind = 'an API key of an OpenAI project'
    model_example = 'gpt-6.1-sol'
    billing = 'openai'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        base_url = (credentials.get('baseUrl') or '').strip()
        m = BASE_URL.fullmatch(base_url) if base_url else None
        if base_url and not m:
            raise ProviderError('The baseUrl container credential is not an OpenAI API address: set '
                                f'https://<region>.api.openai.com, region one of {", ".join(REGIONS)}, for a project '
                                'with data residency, or leave it empty')
        return {'codex_responses': Route(self.endpoint(f'https://{m[1] or ""}api.openai.com/v1' if m else API))}

    def check_key(self, key: str) -> None:
        # A ChatGPT or Codex login and a workload identity federation token are JWTs.
        if key.startswith('eyJ'):
            raise ProviderError('The apiKey container credential is a token, such as a ChatGPT or Codex login: '
                                'set an API key of an OpenAI project')
        if key.startswith('sk-admin-'):
            raise ProviderError('The apiKey container credential is an OpenAI Admin API key, which works only on the '
                                'administration endpoints: set an API key of a project')

    def error_hint(self, status: int | None) -> str:
        return (' (check the apiKey container credential, an API key of an OpenAI project, and for a project with '
                'data residency the baseUrl one)') if status == 401 else ''
