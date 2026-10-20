import re
from typing import ClassVar

from ..hermes import Connection
from .base import Credentials, Provider, ProviderError

# The AWS commercial regions: GovCloud has no global inference profiles, and China has a domain of its own.
REGION = re.compile(r'(af|ap|ca|eu|il|me|mx|sa|us)-[a-z]+-\d{1,2}')
PATHS ={'anthropic_messages': '/anthropic', 'codex_responses': '/openai/v1'}


class BedrockProvider(Provider):
    """Claude and GPT on Amazon Bedrock's bedrock-runtime endpoint with an Amazon Bedrock API key, which Hermes sends
    as x-api-key to the Messages API and as Bearer to the Responses API."""

    name = 'Bedrock'
    # Global cross-Region inference profiles: every commercial region reaches them, at Anthropic's and OpenAI's own
    # prices; a geographic profile costs 10% more. Each model's card in
    # https://docs.aws.amazon.com/bedrock/latest/userguide/models-endpoint-availability.html
    model_ids: ClassVar[dict[str, str]] = {
        **{f'global.anthropic.{m}': m for m in ('claude-haiku-5-5', 'claude-sonnet-5-5', 'claude-opus-5-5',
                                                 'claude-fable-5-1')},
        **{f'global.openai.{m}': m for m in ('gpt-6-astra', 'gpt-6.1-sol', 'gpt-6-sol', 'gpt-6-luna', 'gpt-5.6-sol',
                                              'gpt-5.6-terra', 'gpt-5.6-luna')}}
    default_model, fast_model = 'global.anthropic.claude-sonnet-5-5', 'global.anthropic.claude-haiku-5-5'

    def __init__(self, credentials: Credentials) -> None:
        super().__init__(credentials)
        self._key = self.api_key(credentials, 'an Amazon Bedrock API key')
        region = credentials.get('region') or ''
        if not REGION.fullmatch(region):
            raise ProviderError('The region container credential must name the AWS commercial region the Bedrock API '
                                'key was created in, such as us-east-1')
        host = f'https://bedrock-runtime.{region}.amazonaws.com'
        self._base_urls = {mode: self.endpoint(host + path) for mode, path in PATHS.items()}

    def connection(self, model: str) -> Connection:
        mode = self.api_mode(model)
        return {'provider': 'custom', 'api_mode': mode, 'base_url': self._base_urls[mode], 'api_key': self._key}

    def error_hint(self, status: int | None) -> str:
        return (' (check the apiKey and region container credentials: a Bedrock API key works only in the region it '
                'was created in, and only for the models its IAM user may invoke)') if status in (401, 403) else ''
