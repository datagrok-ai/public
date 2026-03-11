# KnimeLink

Datagrok package for integrating with the KNIME Business Hub.

## What It Does

- Connects to KNIME Business Hub (hub.knime.com or self-hosted) via REST API
- Lists all deployments (REST, data-app, schedule, trigger) for teams the user belongs to
- Passes inputs (tables, parameters, files) and executes workflows
- Displays results as Datagrok DataFrames (REST deployments) or provides link to KNIME Hub (data-app deployments)

## Setup

All settings are configured via standard Datagrok package settings: **Manage > Plugins > KnimeLink**.
If settings are not configured, the app shows a warning toast instead of a non-functional tree item.

| Setting             | Default                  | Description                                 |
|---------------------|--------------------------|---------------------------------------------|
| `baseUrl`           | `https://hub.knime.com`  | The package auto-derives the `api.` prefix  |
| `timeoutMinutes`    | `10`                     | Max time to wait for async job completion   |
| `pollingIntervalMs` | `3000`                   | Interval between job status polls           |

**Credentials** (set in package credentials via Manage > Plugins > KnimeLink > Credentials):
- `appPasswordId` — Application Password ID
- `appPasswordValue` — the generated token

To create an Application Password: click your avatar on hub.knime.com > switch to your **personal** account (not the PRO team) > go to your profile > **Settings** > **Application Passwords** > Create.

**Important**: All deployment types (REST, data-app, schedule, trigger) appear in the tree. Only REST service deployments return API-consumable output data. Data-app deployments produce interactive sessions viewable on KNIME Hub. Raw workflows in your Hub space are not directly executable via API. To deploy: upload workflow to Hub space > create a Deployment.

## Architecture

### Execution Strategy

The package uses different execution strategies based on deployment type (determined by ID prefix):

- **REST services** (`rest:{uuid}`) — **synchronous execution** via `POST /deployments/{id}/execution`.
  Input (tables, variables) is sent in the request body as `InputParameters` (JSON) or as `multipart/form-data`
  when file inputs are present. Results are returned inline.
  This is required because `POST /deployments/{id}/jobs` auto-executes the workflow immediately with
  default data, making it impossible to inject input via the async 2-step flow.
- **Data-apps and other types** — **asynchronous execution**: create job via `POST /deployments/{id}/jobs`,
  then poll via `GET /jobs/{uuid}`. Timeout and polling interval are configurable via package settings.

### Hub Client

`KnimeHubClient` (`src/knime-hub-client.ts`) — single implementation of `IKnimeClient` interface.
- API URL: auto-derived as `api.{hostname}` from the user-entered baseUrl (with fallback for edge cases)
- `GET /accounts/identity` — get current user + teams[]
- `GET /deployments/{scope}` — list deployments for a team (scope = URL-encoded `account:team:{uuid}`)
- `GET /deployments/{id}/open-api` — get workflow OpenAPI spec (input parameter names and types)
- `POST /deployments/{id}/execution` — sync execution for REST services; sends JSON body or `multipart/form-data` when files are present
- `POST /deployments/{id}/jobs` — create async job (auto-executes with defaults, no body input)
- `POST /jobs/{uuid}?async=true` — execute existing job with input body (for async flow)
- `GET /jobs/{uuid}` — poll job status; `workflowState` field: `EXECUTION_FINISHED`, `EXECUTION_FAILED`, etc.; response includes `outputValues`, `outputResources`, `nodeMessages`, and `@controls`
- `GET /jobs/{uuid}/output-resources/{resourceId}` — fetch a named output resource (table/file); handles S3 pre-signed URL redirects via `redirect: 'manual'`
- `DELETE /jobs/{uuid}/execution` — cancel running job
- Auth: Application Password as HTTP Basic
- Team discovery: calls `/accounts/identity`, extracts `teams[].id` for deployment scopes

### Key Files

| File | Purpose |
|------|---------|
| `src/package.ts` | App entry, tree browser, exported functions |
| `src/constants.ts` | `KnimeJobState` enum |
| `src/types.ts` | `KnimeDeployment`, `KnimeInputParam`, `KnimeJobStatus` (with `_rawData` cache), etc. |
| `src/credentials.ts` | `_package.getCredentials()` with session caching |
| `src/knime-client.ts` | `IKnimeClient` interface |
| `src/knime-hub-client.ts` | Business Hub implementation (with `api.` URL derivation) |
| `src/knime-client-factory.ts` | Factory creating `KnimeHubClient` from settings |
| `src/data-conversion.ts` | DataFrame <-> KNIME JSON table format |
| `src/workflow-input-form.ts` | Dynamic input form from workflow parameter descriptors |
| `src/utils.ts` | Async job polling with `DG.TaskBarProgressIndicator`, timeout management |
| `css/knime-link.css` | Styles for workflow view, forms, results, file outputs |

### Tree Browser

Discovers teams via `/accounts/identity`, then lists deployments per team via `GET /deployments/{scope}`.
Each item opens an execution view. The header includes a custom breadcrumb navigation UI.

### Data Conversion

- **DataFrame -> KNIME**: `table-spec`/`table-data` format (column specs as `{name: type}` entries, row data as arrays). Datetime values are formatted as `YYYY-MM-DDTHH:mm:ss.SSS` (no timezone — KNIME `localdatetime` rejects `Z` suffix).
- **KNIME -> DataFrame**: Parses `table-spec`/`table-data` via `knimeSpecDataToDataFrame`, or flat JSON arrays via `DG.DataFrame.fromObjects`.
- **Files**: Sent as multipart/form-data when present (not base64-encoded in JSON)
- **Variables**: Key-value pairs (Container Input Variable)

### Result Visualization

Results are extracted from the sync response or job response:

- **`outputValues`** — inline key-value results; objects with `table-data` become DataFrames, scalars become variable rows
- **`outputResources`** — named resources fetched via `GET /jobs/{uuid}/output-resources/{resourceId}`; text responses (CSV/TSV/JSON) are auto-parsed into DataFrames via `tryParseTextAsDataFrame()`; binary content is shown as a download link with extension guessed from `Content-Type`
- **Output resource pre-fetching** — resources are fetched immediately when a job completes (before any job cleanup) to prevent race conditions with S3 pre-signed URL expiry
- **Error handling** — both `nodeMessages` and `errors` arrays from failed executions are extracted and shown via `grok.shell.error`
- **Flat result parsing** — handles responses where `table-spec`/`table-data` appears at the top level rather than wrapped in `outputValues`/`outputResources`
- If no outputs, shows "Workflow completed with no output."

### Workflow Input Discovery

Input parameters are discovered via `GET /deployments/{id}/open-api`, which returns the workflow's OpenAPI spec.
The parser reads `components.schemas.InputParameters.properties` directly (not from paths — the `/execution`
endpoint references `InputParameters` via `$ref`). File inputs are extracted from the `/execution` endpoint's
`multipart/form-data` schema via `parseMultipartFileInputs()`. Each property is mapped to a Datagrok input
type via `ui.input.forInputType()` using a `KnimeParamType` → `DG.InputType` mapping.

**Input form features:**
- Parameters are grouped into accordion sections when groups are present in the spec
- Input descriptions are shown via info icons with hover tooltips
- Table inputs display column schema in their tooltip
- Default values from the spec are pre-populated in the form
- A "Variables (JSON)" text area is always included for passing flow variables

**Important limitations of the OpenAPI spec:**
- Only **Container Input (Table)**, **Container Input (JSON)**, and **Container Input (File)** nodes are exposed.
- **Container Input (Variable)** nodes are NOT included in the spec. Use the "Variables (JSON)" text area.

When no spec is available, a fallback form with a table picker and "Table Parameter Name" field is shown (default: `table-input`).

### Multipart File Upload

When the workflow has file inputs (detected from the OpenAPI spec's multipart schema), the execution
request is sent as `multipart/form-data` instead of JSON. Non-file inputs (tables, variables) are
serialized to JSON and included as a separate `data` part. This applies to both sync (`/execution`)
and async (`/jobs/{uuid}`) execution paths.

### Exported Functions

- `executeKnimeWorkflow(workflowId, inputJson?, inputTable?, tableParamName?)` — programmatic workflow execution

## Build

```bash
npm install
npm run build    # grok api && grok check --soft && webpack
```

## OpenAPI Specs

The `openapi/` folder contains swagger specs downloaded from `api.hub.knime.com/api-doc/`:
- `accounts-service` — user identity, teams, permissions
- `execution-service` — deployments, jobs, execution contexts, executors

Key schemas: `UserAccount` (includes `teams[]` of `BaseTeam`), `Deployment` (types: rest, schedule, data-app, trigger), `WorkflowJob` (states in `workflowState` field, outputs in `outputValues`/`outputResources`, errors in `nodeMessages`/`errors`, session URL in `@controls`), `InputParameters` (map of param name to value — tables use `table-spec`/`table-data` format).

## Patterns Reused

- App + tree browser decorator pattern from BenchlingLink
- `grok.dapi.fetchProxy()` wrapper from CddVaultLink
- Credential caching from BenchlingLink
