import re

from ..formats import ApiMode
from .base import Credentials, Provider, ProviderError, Route

# A workspace's own URL on AWS, Azure or Google Cloud:
# https://docs.databricks.com/aws/en/workspace/workspace-details (and its Azure and Google Cloud versions)
WORKSPACE = re.compile(r'https://([a-z0-9][a-z0-9-]*\.cloud\.databricks\.com|adb-\d+\.\d+\.azuredatabricks\.net|'
                       r'\d+\.\d+\.gcp\.databricks\.com)/?')


class DatabricksProvider(Provider):
    """Claude and GPT through a Databricks workspace's AI Gateway, with a Databricks token as Bearer on both routes
    (https://docs.databricks.com/aws/en/machine-learning/model-serving/query-anthropic-messages,
    …/query-openai-responses). A model is named as the gateway serves it, such as system.ai.claude-sonnet-5-5."""

    name = 'Databricks'
    key_kind = 'a Databricks personal access token'
    model_example = 'system.ai.claude-sonnet-5-5'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        workspace = credentials.get('baseUrl') or ''
        if not WORKSPACE.fullmatch(workspace):
            raise ProviderError('The baseUrl container credential is not a Databricks workspace URL: set '
                                'https://<workspace>.cloud.databricks.com, https://adb-<id>.<n>.azuredatabricks.net '
                                'or https://<id>.<n>.gcp.databricks.com')
        base = workspace.rstrip('/')
        return {'anthropic_messages': Route(self.endpoint(base + '/ai-gateway/anthropic'), bearer=True),
                'codex_responses': Route(self.endpoint(base + '/ai-gateway/openai/v1'), bearer=True)}

    def error_hint(self, status: int | None) -> str:
        return (' (check the apiKey container credential: a token of this workspace whose user or service principal '
                'may query the model)') if status in (401, 403) else ''
