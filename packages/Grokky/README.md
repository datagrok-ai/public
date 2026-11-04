# Grokky

`Grokky` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform

## Model credentials

The AI model Grokky runs on is set in the package's credentials: **Browse > Platform > Plugins**, right-click
**Grokky**, **Credentials...**. **All users** holds the default; pick a group on the left to give its members
their own key and model. Keep **Allow client side** off, so the values reach only the platform and Grokky's
runtime container, never the browser.

| Key         | Meaning                                                                                        |
|-------------|------------------------------------------------------------------------------------------------|
| `provider`  | `Anthropic` (default), `OpenAI`, `OpenRouter`, `Azure`, `Bedrock`, `Databricks` or `Custom`    |
| `apiKey`    | The provider's key; required                                                                   |
| `model`     | The model every turn runs on, as the provider names it (e.g. `claude-sonnet-5-5`); required    |
| `baseUrl`   | Azure: the resource; Databricks: the workspace; Custom: the endpoint; OpenAI: a data-residency host |
| `region`    | Bedrock: the AWS region the key was created in                                                 |
| `apiFormat` | Custom: `anthropic-messages`, `openai-responses` or `openai-chat-completions`                  |
| `service`   | Databricks: a model service of the workspace to call `model` through                           |

A user's chat runs on the credentials of their group; a user the package has none for gets an error. Leave the
credentials of the `grokky-claude-runtime` container itself empty: the platform lays them over the package's.
A change reaches a user when their AI panel next connects (a page reload). `POST /credentials/for/Grokky?group=<group>` saves a group's values with **Allow client side** on:
turn it off in the dialog afterwards.
