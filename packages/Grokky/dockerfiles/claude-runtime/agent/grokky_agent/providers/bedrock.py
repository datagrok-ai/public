import re

from ..formats import ApiMode
from ..hermes import Connection
from .base import CLAUDE, Credentials, Provider, ProviderError, Route

# The AWS commercial regions: GovCloud has no global inference profiles, and China has a domain of its own.
REGION = re.compile(r'(af|ap|ca|eu|il|me|mx|sa|us)-[a-z]+-\d{1,2}')


class BedrockProvider(Provider):
    """Models on Amazon Bedrock's bedrock-runtime endpoint with an Amazon Bedrock API key: Claude on the Messages API
    (the key as x-api-key), GPT on the Responses API (as Bearer), and every other model (Qwen, Kimi, ...) on the
    Converse API. A model is named by its inference profile,
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
                'codex_responses': Route(self.endpoint(host + '/openai/v1')),
                'bedrock_converse': Route(self.endpoint(host))}

    def mode(self) -> ApiMode:
        if CLAUDE.search(self.model):
            return 'anthropic_messages'
        return 'codex_responses' if re.match(r'([a-z]+\.)?openai\.', self.model) else 'bedrock_converse'

    def connection(self) -> Connection:
        if self.format.mode != 'bedrock_converse':
            return super().connection()
        # Hermes' own Bedrock provider; hermes.use_bedrock_key gives its boto3 clients the key.
        return {'provider': 'bedrock', 'api_mode': 'bedrock_converse', 'base_url': self.route.url, 'api_key': 'aws-sdk'}

    def error_hint(self, status: int | None) -> str:
        return (' (check the apiKey and region container credentials: a Bedrock API key works only in the region it '
                'was created in, and only for the models its IAM user may invoke)') if status in (401, 403) else ''
