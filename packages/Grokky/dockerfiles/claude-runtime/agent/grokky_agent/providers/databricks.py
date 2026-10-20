import re
from typing import ClassVar

from ..formats import ApiMode
from ..hermes import Connection
from .base import Credentials, Provider, ProviderError

# A workspace's own URL on AWS, Azure or Google Cloud:
# https://docs.databricks.com/aws/en/workspace/workspace-details (and its Azure and Google Cloud versions)
WORKSPACE = re.compile(r'https://([a-z0-9][a-z0-9-]*\.cloud\.databricks\.com|adb-\d+\.\d+\.azuredatabricks\.net|'
                       r'\d+\.\d+\.gcp\.databricks\.com)/?')
# https://docs.databricks.com/aws/en/machine-learning/model-serving/query-anthropic-messages, …/query-openai-responses
ROUTES: dict[ApiMode, str] = {'anthropic_messages': '/ai-gateway/anthropic', 'codex_responses': '/ai-gateway/openai/v1'}
# Databricks hosts Claude Sonnet 5.5 for text input only, under its model service and its pay-per-token endpoint name.
TEXT_ONLY = {'system.ai.claude-sonnet-5-5', 'databricks-claude-sonnet-5-5'}


class DatabricksProvider(Provider):
    """Claude and GPT through a Databricks workspace's AI Gateway, with a Databricks token."""

    name = 'Databricks'
    # Unity Gateway's model service names, which write the version's dots as hyphens:
    # https://docs.databricks.com/aws/en/machine-learning/foundation-model-apis/supported-models
    model_ids: ClassVar[dict[str, str]] = {f'system.ai.{m.replace(".", "-")}': m for m in (
        'claude-haiku-5-5', 'claude-sonnet-5-5', 'claude-opus-5-5', 'claude-fable-5-1', 'gpt-6-luna', 'gpt-6-sol',
        'gpt-6.1-sol', 'gpt-6-astra', 'gpt-5.6-luna', 'gpt-5.6-terra', 'gpt-5.6-sol')}
    default_model, fast_model = 'system.ai.claude-sonnet-5-5', 'system.ai.claude-haiku-5-5'

    def __init__(self, credentials: Credentials) -> None:
        super().__init__(credentials)
        workspace = credentials.get('baseUrl') or ''
        if not WORKSPACE.fullmatch(workspace):
            raise ProviderError('The baseUrl container credential is not a Databricks workspace URL: set '
                                'https://<workspace>.cloud.databricks.com, https://adb-<id>.<n>.azuredatabricks.net '
                                'or https://<id>.<n>.gcp.databricks.com')
        self._urls = {mode: self.endpoint(workspace.rstrip('/') + path) for mode, path in ROUTES.items()}
        key = self.api_key(credentials, 'a Databricks personal access token')
        # Hermes sends a callable's key as Bearer, which both routes take; a string goes to Claude's as x-api-key.
        self._key = lambda: key

    def connection(self, model: str) -> Connection:
        api_mode = self.api_mode(model)
        return {'provider': 'custom', 'api_mode': api_mode, 'base_url': self._urls[api_mode], 'api_key': self._key}

    def images(self, model: str) -> bool:
        return model not in TEXT_ONLY and super().images(model)

    def error_hint(self, status: int | None) -> str:
        return (' (check the apiKey container credential: a token of this workspace whose user or service principal '
                'may query the model)') if status in (401, 403) else ''
