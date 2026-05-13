---
name: rest-api
description: Call the Datagrok REST API from outside the platform — authenticate, transfer files/tables, build dashboards, invoke functions
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# rest-api

## When to use

You are integrating with a Datagrok server from a non-Datagrok process —
a CI pipeline, a Python notebook, a Node service, a curl-based smoke test
— and need to upload tables, fetch files, build dashboards, or invoke
server-side functions. Triggers: "call Datagrok from a script",
"automate dashboard creation", "POST a CSV from a cron job", "OpenAPI
spec". For *in-package* HTTP, use `grok.dapi.fetchProxy` (see
`access-data` skill) — this skill is for external callers.

## Prerequisites

- A Datagrok deployment URL (e.g. `https://public.datagrok.ai`). The
 REST base URL is that URL plus `/api` (knowledge `DG-FACT-230`).
- A personal API key from the user profile page (`<GROK_HOST>/u`,
 "API Key" link). Every request sends it as
 `Authorization: Bearer <KEY>` — caller adds the `Bearer ` prefix
 literally, including the trailing space (knowledge `DG-FACT-229`).
- `curl`, or the official Python client
 (`pip install datagrok-api`) — same surface, less boilerplate.

## Steps

1. **Set the base URL and API key in the shell.**
 ```bash
 export GROK_HOST=https://public.datagrok.ai
 export GROK_API=$GROK_HOST/api
 export API_KEY="Bearer xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
 ```
 Expected: every later step substitutes `$GROK_API` and `$API_KEY`. The
 `Bearer ` prefix lives inside `$API_KEY` so headers stay one-liners
 (knowledge `DG-FACT-229`).

2. **Smoke-test auth by fetching the OpenAPI spec.**
 The spec is the authoritative endpoint inventory; the help article
 lists only a curated subset (knowledge `DG-FACT-234`).
 ```bash
 curl -fsS -H "Authorization: $API_KEY" "$GROK_API/public/api.yaml" \
 | head -20
 ```
 Expected: YAML starting with `openapi: 3...`. A `401` means the key
 or `Bearer ` prefix is wrong; a `404` means `$GROK_API` is missing
 `/api`.

3. **Download a table as CSV.**
 Tables are addressed by GUID *or* Grok name — but Grok names use `:`
 as separator and the REST API requires `.` instead. Replace before
 sending (knowledge `DG-FACT-232`).
 ```bash
 TABLE="JohnDoe.MyProject.Cars" # NOT JohnDoe:MyProject:Cars
 curl -fsS -H "Authorization: $API_KEY" \
 "$GROK_API/public/v1/tables/$TABLE" -o cars.csv
 ```
 Expected: `cars.csv` on disk; `head cars.csv` shows the column
 header. Endpoint shape: `GET /public/v1/tables/{name}` (knowledge
 `DG-FACT-231`, `DG-FACT-233`).

4. **Upload a CSV as a table.**
 `POST` body is the CSV bytes; `Content-Type: text/csv`. If the path
 includes a project segment, the table is added to that project; if
 the table name already exists, its data is replaced.
 ```bash
 curl -fsS -X POST \
 -H "Authorization: $API_KEY" \
 -H "Content-Type: text/csv" \
 --data-binary @iris.csv \
 "$GROK_API/public/v1/tables/JohnDoe.iris"
 ```
 Expected: JSON `{"ID": "<guid>",...}` echoing identifiers for the
 newly uploaded table.

5. **Upload/download a file under a connector.**
 File paths are `<connector>/<path>`. The connector identifier obeys
 the same colon→period rule (knowledge `DG-FACT-232`). Body is binary;
 `Content-Type: application/octet-stream`.
 ```bash
 CONN="JohnDoe.Home" # was JohnDoe:Home
 curl -fsS -X POST \
 -H "Authorization: $API_KEY" \
 -H "Content-Type: application/octet-stream" \
 --data-binary @report.pdf \
 "$GROK_API/public/v1/files/$CONN/reports/report.pdf"
 ```
 Expected: HTTP 200; the file appears under
 `<GROK_HOST>/u/files/$CONN/reports/`.

6. **Create a dashboard from uploaded tables.**
 `POST /public/v1/dashboards/{name}/{table_ids}`. `table_ids` is
 comma-separated; the optional JSON body is a project layout
 (knowledge `DG-FACT-233`).
 ```bash
 curl -fsS -X POST \
 -H "Authorization: $API_KEY" \
 -H "Content-Type: application/json" \
 -d @layout.json \
 "$GROK_API/public/v1/dashboards/MyDashboard/$TABLE_ID_1,$TABLE_ID_2"
 ```
 Expected: JSON with the new dashboard's id; visit
 `<GROK_HOST>/p/<dashboard-id>` to confirm.

7. **Invoke a server-side function.**
 The article documents `POST /public/v1/{name}/call`, but that path
 404s — the actual segment is `/public/v1/functions/{name}/call`. Body is the parameter object as JSON.
 ```bash
 curl -fsS -X POST \
 -H "Authorization: $API_KEY" \
 -H "Content-Type: application/json" \
 -d '{"a": 1, "b": 2}' \
 "$GROK_API/public/v1/functions/JohnDoe.MyFunction/call"
 ```
 Expected: JSON-encoded function return value.

8. **Skip the curl boilerplate with the Python client.**
 `pip install datagrok-api` then construct a `DatagrokClient` once;
 resources cover files, tables, dashboards, functions, connections,
 users, groups, shares (knowledge `DG-FACT-231`).
 ```python
 from datagrok_api import DatagrokClient
 grok = DatagrokClient(base_url="https://public.datagrok.ai/api",
 api_key="Bearer xxxx")
 df = grok.tables.download("JohnDoe:MyTable") # ':' OK — client substitutes
 grok.functions.call("JohnDoe:MyFunction", {"a": 1, "b": 2})
 ```
 Expected: the client replaces `:` with `.` internally
 (knowledge `DG-FACT-232`) and prefixes `Bearer ` is the caller's job
 (knowledge `DG-FACT-229`).

## Common failure modes

- **`401 Unauthorized` on every call.** The API key is missing the
 literal `Bearer ` prefix (with a trailing space) or `$API_KEY` is
 unquoted and the shell ate the space (knowledge `DG-FACT-229`). Fix:
 `API_KEY="Bearer xxxx"` and always quote `"$API_KEY"` in headers.
- **`404 Not Found` calling a function.** You followed the article path
 `/public/v1/{name}/call` literally. Fix:
 insert `functions/` — `/public/v1/functions/{name}/call`.
- **`404 Not Found` on tables/files with the right name.** The
 identifier still has `:` separators; REST API requires `.`
 (knowledge `DG-FACT-232`). Fix:
 `name=$(echo "$NAME" | tr ':' '.')` before composing the URL.
- **Wrong base URL — endpoints route to the web UI.** `$GROK_HOST`
 works in the browser but the REST root is `$GROK_HOST/api`
 (knowledge `DG-FACT-230`). Fix: define `GROK_API="$GROK_HOST/api"`
 and use it for every REST call.
- **Successful HTTP 200 but the response body says `ApiError...`.**
 Datagrok signals application-level errors with an `api-error` header
 even on 2xx. Fix: check `response.headers['api-error']` and surface
 the body — the Python client and the Grokky package both do this
 (`python-api/datagrok_api/http_client.py:48`,
 `packages/Grokky/dockerfiles/claude-runtime/src/shared-api-client.ts:51`).
- **CSV upload arrives mangled.** `curl` defaulted to `--data` which
 strips newlines and applies form encoding. Fix: use `--data-binary
 @file.csv` plus `Content-Type: text/csv`.

## See also

- Source articles:
 - `help/develop/packages/rest-api.md`
 - `help/develop/server-management.md` (CLI for the same surface)
- Knowledge:
 - `docs/_internal/knowledge/knowledge-graph.md` — facts
 `DG-FACT-229` (auth header), `DG-FACT-230` (base URL),
 `DG-FACT-231` (`/public/v1/` prefix), `DG-FACT-232` (colons →
 periods), `DG-FACT-233` (entity surfaces), `DG-FACT-234` (OpenAPI
 URL), and drift (function-call path).
- Related skills:
 - `access-data` — for *in-package* HTTP via `grok.dapi.fetchProxy`
 (CORS bypass + caching) when the caller already runs inside
 Datagrok.
 - `js-api` — `grok.functions.call` is the in-platform equivalent of
 step 7's REST endpoint.
