import re

from ..formats import ApiMode
from .base import Credentials, Provider, ProviderError, Route

# The AWS commercial regions: GovCloud has no global inference profiles, and China has a domain of its own.
REGION = re.compile(r'(af|ap|ca|eu|il|me|mx|sa|us)-[a-z]+-\d{1,2}')


class BedrockProvider(Provider):
    """Claude and GPT on Amazon Bedrock's bedrock-runtime endpoint with an Amazon Bedrock API key, which Hermes sends
    as x-api-key to the Messages API and as Bearer to the Responses API. A model is named by its inference profile,
    such as global.anthropic.claude-sonnet-5-5: a global profile reaches the model from every commercial region at the
    vendor's own price, a geographic one costs 10% more
    (https://docs.aws.amazon.com/bedrock/latest/userguide/models-endpoint-availability.html)."""

    name = 'Bedrock'
    key_kind = 'an Amazon Bedrock API key'
    model_example = 'global.anthropic.claude-sonnet-5-5'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        region = credentials.get('region') or ''
        if not REGION.fullmatch(region):
            raise ProviderError('The region container credential must name the AWS commercial region the Bedrock API '
                                'key was created in, such as us-east-1')
        host = f'https://bedrock-runtime.{region}.amazonaws.com'
        return {'anthropic_messages': Route(self.endpoint(host + '/anthropic')),
                'codex_responses': Route(self.endpoint(host + '/openai/v1'))}

    def error_hint(self, status: int | None) -> str:
        return (' (check the apiKey and region container credentials: a Bedrock API key works only in the region it '
                'was created in, and only for the models its IAM user may invoke)') if status in (401, 403) else ''
