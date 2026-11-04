# Grokky

`Grokky` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform

## Model credentials

Grokky's model is set in the package's credentials: **Browse > Platform > Plugins**, right-click **Grokky**,
**Credentials...**. **All users** is the default; a group on the left gives its members their own key and model.
Keep **Allow client side** off, so the values never reach the browser.

| Key         | Meaning                                                                                        |
|-------------|------------------------------------------------------------------------------------------------|
| `provider`  | `Anthropic` (default), `OpenAI`, `OpenRouter`, `Azure`, `Bedrock`, `Databricks` or `Custom`    |
| `apiKey`    | The provider's key; required                                                                   |
| `model`     | The model every turn runs on, as the provider names it (e.g. `claude-sonnet-5-5`); required    |
| `baseUrl`   | Azure: the resource; Databricks: the workspace; Custom: the endpoint; OpenAI: a data-residency host |
| `region`    | Bedrock: the AWS region the key was created in; required for Bedrock                           |
| `apiFormat` | Custom: `anthropic-messages`, `openai-responses` or `openai-chat-completions`                  |
| `service`   | Databricks: a model service of the workspace (`<catalog>.<schema>.<name>`) to call `model` through |

- A user without credentials gets an error in the AI panel.
- Changes reach a user on their next connection (a page reload).
- Do not set provider credentials on the `grokky-claude-runtime` container: the platform lays them over the package's.
- `POST /credentials/for/Grokky?group=<group>` saves with **Allow client side** on: turn it off in the dialog.
