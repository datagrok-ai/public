import re

from ..formats import ApiMode
from ..protocol import JsonObject
from .base import Credentials, Provider, ProviderError, Route

# A workspace's own URL on AWS, Azure or Google Cloud:
# https://docs.databricks.com/aws/en/workspace/workspace-details (and its Azure and Google Cloud versions)
WORKSPACE = re.compile(r'https://([a-z0-9][a-z0-9-]*\.cloud\.databricks\.com|adb-\d+\.\d+\.azuredatabricks\.net|'
                       r'\d+\.\d+\.gcp\.databricks\.com)/?')
# A Unity Catalog model service: <catalog>.<schema>.<name>.
SERVICE = re.compile(r'[\w-]+\.[\w-]+\.[\w-]+')


class DatabricksProvider(Provider):
    """Claude and GPT through a Databricks workspace's AI Gateway, with a Databricks token as Bearer on both routes
    (https://docs.databricks.com/aws/en/machine-learning/model-serving/query-anthropic-messages,
    …/query-openai-responses). A model is named as the gateway serves it, such as system.ai.claude-sonnet-5-5; with
    the `service` credential, a model service of the workspace's own (an external provider's model) is named in the
    Databricks-Model-Provider-Service header, and `model` is the provider's model id the service allows."""

    name = 'Databricks'
    key_kind = 'a Databricks personal access token'
    model_example = 'system.ai.claude-sonnet-5-5'

    def routes(self, credentials: Credentials) -> dict[ApiMode, Route]:
        workspace = credentials.get('baseUrl') or ''
        if not WORKSPACE.fullmatch(workspace):
            raise ProviderError('The baseUrl container credential is not a Databricks workspace URL: set '
                                'https://<workspace>.cloud.databricks.com, https://adb-<id>.<n>.azuredatabricks.net '
                                'or https://<id>.<n>.gcp.databricks.com')
        service = (credentials.get('service') or '').strip()
        if service and not SERVICE.fullmatch(service):
            raise ProviderError('The service container credential must name a model service as '
                                '<catalog>.<schema>.<name>, such as workspace.default.claude-haiku-5-5-ext')
        self.service = service
        base = workspace.rstrip('/')
        return {'anthropic_messages': Route(self.endpoint(base + '/ai-gateway/anthropic'), bearer=True),
                'codex_responses': Route(self.endpoint(base + '/ai-gateway/openai/v1'), bearer=True)}

    def patch_request(self, api_mode: str, request: JsonObject, system_prompt: str,
                      output_schema: JsonObject | None) -> JsonObject:
        out = super().patch_request(api_mode, request, system_prompt, output_schema)
        if self.service:
            out['extra_headers'] = {**(out.get('extra_headers') or {}),
                                    'Databricks-Model-Provider-Service': self.service}
        return out

    def error_hint(self, status: int | None) -> str:
        if status in (401, 403):
            return (' (check the apiKey container credential: a token of this workspace whose user or service '
                    'principal may query the model' +
                    (', and that model is a model id the service allows)' if self.service else ')'))
        return ' (check the model and service container credentials)' if status == 404 and self.service else ''
