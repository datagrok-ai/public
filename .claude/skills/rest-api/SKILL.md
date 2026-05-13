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

1. **Set the base URL and API key in the shell.** `Bearer ` prefix lives inside `$API_KEY` (`DG-FACT-229`).
 ```bash
 export GROK_HOST=https://public.datagrok.ai
 export GROK_API=$GROK_HOST/api
 export API_KEY="Bearer xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
 ```

2. **Smoke-test auth by fetching the OpenAPI spec** — the authoritative endpoint inventory (`DG-FACT-234`).
 ```bash
 curl -fsS -H "Authorization: $API_KEY" "$GROK_API/public/api.yaml" \
 | head -20
 ```
 `401` = key/prefix wrong; `404` = `$GROK_API` missing `/api`.

3. **Download a table as CSV.** Tables are addressed by GUID or Grok name; Grok names use `:` but REST requires `.` (`DG-FACT-232`).
 ```bash
 TABLE="JohnDoe.MyProject.Cars" # NOT JohnDoe:MyProject:Cars
 curl -fsS -H "Authorization: $API_KEY" \
 "$GROK_API/public/v1/tables/$TABLE" -o cars.csv
 ```
 Endpoint: `GET /public/v1/tables/{name}` (`DG-FACT-231`, `DG-FACT-233`).

4. **Upload a CSV as a table.** Body is raw CSV with `Content-Type: text/csv`; existing names are replaced.
 ```bash
 curl -fsS -X POST \
 -H "Authorization: $API_KEY" \
 -H "Content-Type: text/csv" \
 --data-binary @iris.csv \
 "$GROK_API/public/v1/tables/JohnDoe.iris"
 ```
 Expected: JSON `{"ID": "<guid>",...}` echoing identifiers for the
 newly uploaded table.

5. **Upload/download a file under a connector.** Path is `<connector>/<path>`; same colon→period rule (`DG-FACT-232`). Body is binary with `Content-Type: application/octet-stream`.
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

6. **Create a dashboard from uploaded tables** — `POST /public/v1/dashboards/{name}/{table_ids}` (comma-separated; optional JSON-body layout). See `DG-FACT-233`.
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

8. **Skip the curl boilerplate with the Python client.** `pip install datagrok-api`; resources cover files/tables/dashboards/functions/connections/users/groups/shares (`DG-FACT-231`). Client substitutes `:`→`.` (`DG-FACT-232`); caller still adds `Bearer ` (`DG-FACT-229`).
 ```python
 from datagrok_api import DatagrokClient
 grok = DatagrokClient(base_url="https://public.datagrok.ai/api",
 api_key="Bearer xxxx")
 df = grok.tables.download("JohnDoe:MyTable") # ':' OK — client substitutes
 grok.functions.call("JohnDoe:MyFunction", {"a": 1, "b": 2})
 ```

## Common failure modes

- **`401 Unauthorized`** — missing `Bearer ` prefix, or unquoted `$API_KEY` (`DG-FACT-229`). Always quote `"$API_KEY"`.
- **`404` calling a function** — article path is wrong; insert `functions/`: `/public/v1/functions/{name}/call`.
- **`404` on tables/files with the right name** — identifier still has `:`; replace with `.` (`DG-FACT-232`).
- **Wrong base URL routes to web UI** — REST root is `$GROK_HOST/api` (`DG-FACT-230`).
- **HTTP 200 but body says `ApiError...`** — check `response.headers['api-error']` (errors surface even on 2xx).
- **CSV upload mangled** — use `--data-binary @file.csv` not `--data`.

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
