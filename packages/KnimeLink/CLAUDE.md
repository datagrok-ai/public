# KnimeLink

Datagrok package for integrating with the KNIME Business Hub.

## What It Does

- Connects to KNIME Business Hub (hub.knime.com or self-hosted) via REST API
- Lists REST service deployments for teams the user belongs to
- Passes inputs (tables, parameters, files) and executes workflows
- Displays results as Datagrok DataFrames (REST deployments) or provides link to KNIME Hub (data-app deployments)

## Maintenance Rule

When making significant logic changes (new files, changed architecture, renamed conventions,
altered data flow, new/removed functions), update this CLAUDE.md to reflect the change before
considering the task complete.

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

**Important**: Only REST service deployments appear in the tree and are executable via the API. Raw workflows in your Hub space are not directly executable via API. To deploy: upload workflow to Hub space > create a REST Deployment.

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
- `GET /knime/rest/v4/repository/{workflowId}:image` — workflow SVG diagram (catalog service); loaded via `<img>` element, not `fetch()`, because the endpoint redirects to S3 and all `fetch()` approaches fail (see `getWorkflowImageUrl()` JSDoc for details)
- `POST /deployments/{id}/execution` — sync execution for REST services; sends JSON body or `multipart/form-data` when files are present
- `POST /deployments/{id}/jobs` — create async job (auto-executes with defaults, no body input)
- `POST /jobs/{uuid}?async=true` — execute existing job with input body (for async flow)
- `GET /jobs/{uuid}` — poll job status; `workflowState` field: `EXECUTION_FINISHED`, `EXECUTION_FAILED`, etc.; response includes `outputValues`, `outputResources`, `nodeMessages`, and `@controls`
- `GET /jobs/{uuid}/output-resources/{resourceId}` — fetch a named output resource (table/file); handles S3 pre-signed URL redirects via URL-embedded credentials + `redirect: 'follow'`
- `DELETE /jobs/{uuid}/execution` — cancel running job
- Auth: Application Password as HTTP Basic
- Team discovery: calls `/accounts/identity`, extracts `teams[].id` for deployment scopes

### Key Files

| File | Purpose |
|------|---------|
| `src/package.ts` | App entry, tree browser, exported functions |
| `src/constants.ts` | `KnimeJobState` enum |
| `src/types.ts` | `KnimeDeployment` (includes `workflowId`), `KnimeInputParam`, `KnimeJobStatus` (with `_rawData` cache), etc. |
| `src/credentials.ts` | `_package.getCredentials()` with session caching |
| `src/knime-client.ts` | `IKnimeClient` interface |
| `src/knime-hub-client.ts` | Business Hub implementation (with `api.` URL derivation) |
| `src/knime-client-factory.ts` | Factory creating `KnimeHubClient` from settings |
| `src/function-registry.ts` | Dynamic Datagrok function registration from KNIME deployment specs |
| `src/function-cache.ts` | Persistent cache via `grok.userSettings`: sync load at startup, async background refresh with diff |
| `src/data-conversion.ts` | DataFrame <-> KNIME JSON table format |
| `src/utils.ts` | Async job polling with `DG.TaskBarProgressIndicator`, timeout management |
| `css/knime-link.css` | Styles for workflow preview, breadcrumbs, progress indicators, results display |

### Tree Browser

Uses a two-phase approach:
- **Phase 1 (cached)**: Loads cached entries via `loadCachedEntries()` and creates tree items immediately
  with already-registered functions (from autostart). Items are usable before the network request completes.
- **Phase 2 (live refresh)**: Shows a `TaskBarProgressIndicator`, fetches live deployments from the API,
  updates existing items' values with fresh functions, adds new items, and removes items whose deployment
  names are no longer found on the Hub.

Test deployments (prefixed with `datagrok_test_`) are always grouped last under a "Test workflows" node.

Selecting a deployment registers it as a Datagrok function and shows two previews:
- **Function preview** (`grok.shell.preview`) — the platform's built-in function preview with inputs/outputs
- **Workflow image** (`grok.shell.o`, context panel) — SVG diagram loaded via an `<img>` element using
  `getWorkflowImageUrl()`, which builds a Datagrok proxy URL with KNIME credentials embedded in the URL.
  The `<img>` approach is required because the KNIME catalog `:image` endpoint returns a 307 redirect to an
  S3 pre-signed URL, and all `fetch()`-based approaches fail (see `getWorkflowImageUrl()` JSDoc).

The `workflowId` (repository item ID) is captured from the deployment metadata.
Test deployments (prefixed with `datagrok_test_`) are grouped under a "Test workflows" node.
The header includes a custom breadcrumb navigation UI.

### Data Conversion

- **DataFrame -> KNIME**: `table-spec`/`table-data` format (column specs as `{name: type}` entries, row data as arrays). Datetime values are formatted as `YYYY-MM-DDTHH:mm:ss.SSS` (no timezone — KNIME `localdatetime` rejects `Z` suffix).
- **KNIME -> DataFrame**: Parses `table-spec`/`table-data` via `knimeSpecDataToDataFrame`, or flat JSON arrays via `DG.DataFrame.fromObjects`.
- **Files**: Sent as multipart/form-data when present (not base64-encoded in JSON)
- **Variables**: Key-value pairs (Container Input Variable)

### Result Visualization

Results are extracted from the sync response or job response:

- **`outputValues`** — inline key-value results; objects with `table-data` become DataFrames, scalars become variable rows
- **`outputResources`** — named resources fetched via `GET /jobs/{uuid}/output-resources/{resourceId}`; plain string values are returned as scalars; JSON arrays are parsed via `knimeTableToDataFrame()`, CSV/TSV text is parsed via `DG.DataFrame.fromCsv()`; binary content is returned as a blob
- **Output resource pre-fetching** — resources are fetched immediately when a job completes (before any job cleanup) to prevent race conditions with S3 pre-signed URL expiry
- **Error handling** — both `nodeMessages` and `errors` arrays from failed executions are extracted; errors throw exceptions
- **Flat result parsing** — handles responses where `table-spec`/`table-data` appears at the top level rather than wrapped in `outputValues`/`outputResources`

### Workflow Input Discovery

Input parameters are discovered via `GET /deployments/{id}/open-api`, which returns the workflow's OpenAPI spec.
The parser reads `components.schemas.InputParameters.properties` directly (not from paths — the `/execution`
endpoint references `InputParameters` via `$ref`). File inputs are extracted from the `/execution` endpoint's
`multipart/form-data` schema via the private `parseMultipartFileInputs()` method. Each property is mapped
to a `KnimeParamType` and then to a Datagrok function parameter type via the `knimeInputTypeToDgType` mapping.
There is no custom input form — inputs are exposed as Datagrok function parameters via `grok.functions.register()`,
and the platform's built-in function UI handles the form rendering.

**Important limitations of the OpenAPI spec:**
- Only **Container Input (Table)**, **Container Input (JSON)**, and **Container Input (File)** nodes are exposed.
- **Container Input (Variable)** nodes are NOT included in the spec.

When no spec is available (fetch fails), the function is registered with no input parameters.

### Multipart File Upload

When the workflow has file inputs (detected from the OpenAPI spec's multipart schema), the execution
request is sent as `multipart/form-data` instead of JSON. Non-file inputs (tables, variables) are
serialized to JSON and included as a separate `data` part. This applies to both sync (`/execution`)
and async (`/jobs/{uuid}`) execution paths.

### Dynamic Function Registration

`function-registry.ts` dynamically registers Datagrok functions from KNIME deployment OpenAPI specs at runtime
via `grok.functions.register()`. When a deployment is selected in the tree browser, `getOrRegisterFunc()`:

1. Fetches the workflow's OpenAPI spec via `client.getWorkflowInputs()`
2. Builds `ParamMeta[]` from the parsed `KnimeInputParam[]`, sanitizing names for Datagrok compatibility
3. Builds `OutputMeta[]` from parsed `KnimeOutputParam[]`, with sanitized names and Datagrok types
4. Generates a function signature string (e.g., `object MyWorkflow(dataframe input_table, string threshold)`)
5. Creates a run callback that converts Datagrok types to KNIME formats, executes, and extracts results
6. Registers the function with `grok.functions.register()`, passing `outputs` array for multi-output workflows, and sets descriptions/defaults on input properties

**Multiple outputs**: Workflows with multiple outputs (tables, scalars, resources) register each output as a
separate named output parameter via the `outputs` field of `grok.functions.register()`. The run callback
returns a named object `{outputName: value}` keyed by sanitized output names. Single-output workflows
return the value directly for backward compatibility.

**Grouped parameters**: When a KNIME input is an object with typed sub-properties, each sub-property becomes
a separate function input with a sanitized name `{group}_{param}`. The description for grouped inputs
includes `Part of "{group}" parameter` (with the parent's description in parentheses if available); the input's own description is prepended before the group label, so users
understand which inputs belong together. At execution time, grouped inputs are reassembled into a nested
object under the original parent key before being sent to KNIME.

**Input/output type mapping**: Input parameters use `knimeInputTypeToDgType` (KNIME `file` → DG `file`),
while output parameters use `knimeOutputTypeToDgType` (KNIME `file` → DG `string`). File outputs map to
`string` because KNIME returns filenames, not `DG.FileInfo` objects, and the platform's `FileInput` widget
would crash trying to call `get$fullPath` on a plain string.

**Name sanitization**: Function names are prefixed with `Knime_` and stripped of non-alphanumeric characters
(e.g., `My Workflow` → `Knime_My_Workflow`). Parameter names are stripped similarly; numeric-leading
param names get a `p_` prefix; duplicate names get a numeric suffix.

### Persistent Function Cache

The autostart function (`knimeLinkAutostart`) registers KNIME workflow functions in two phases:

**Phase 1 — instant from `grok.userSettings` (synchronous, no network)**:
`loadCachedEntries()` reads deployment specs from `grok.userSettings` (storage: `KnimeLinkFuncCache`).
This is a synchronous in-memory read — functions appear instantly on startup with zero async calls.
Each entry stores `KnimeDeployment`, `KnimeWorkflowSpec`, and a pre-computed signature for diff detection.

**Phase 2 — background refresh (async, network)**:
`refreshAndUpdateCache()` fetches live deployments and specs from KNIME Hub, computes signature diffs
against cached entries, registers new/changed functions, and writes updated entries back to `grok.userSettings`
so the next startup's Phase 1 has fresh data. Per-deployment spec fetch failures are handled individually —
the cached version is kept.

### Exported Functions

The package exports three decorated functions:
- `knimeLinkAutostart()` — autostart function that registers KNIME functions from cache, then refreshes in background
- `knimeLinkApp()` — app entry point at Browse > Compute > KNIME
- `knimeLinkAppTreeBrowser()` — tree browser that lists REST deployments per team

## Build

```bash
npm install
npm run build    # grok api && grok check --soft && webpack
```

## OpenAPI Specs

The `openapi/` folder contains swagger specs downloaded from `api.hub.knime.com/api-doc/`:
- `accounts-service` — user identity, teams, permissions
- `execution-service` — deployments, jobs, execution contexts, executors
- `catalog-service` — repository items, workflow images (`GET /repository/{id}:image`)

Key schemas: `UserAccount` (includes `teams[]` of `BaseTeam`), `Deployment` (types: rest, schedule, data-app, trigger), `WorkflowJob` (states in `workflowState` field, outputs in `outputValues`/`outputResources`, errors in `nodeMessages`/`errors`, session URL in `@controls`), `InputParameters` (map of param name to value — tables use `table-spec`/`table-data` format).

## Patterns Reused

- App + tree browser decorator pattern from BenchlingLink
- `grok.dapi.fetchProxy()` wrapper from CddVaultLink
- Credential caching from BenchlingLink
