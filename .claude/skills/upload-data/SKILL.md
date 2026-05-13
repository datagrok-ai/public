---
name: upload-data
version: 0.1.0
description: |
  Push a CSV into a Datagrok workspace from OUTSIDE Datagrok — a shell
  script, cron job, CI pipeline, lab instrument adapter, or standalone
  Python notebook — by POSTing to the legacy `/api/files/<LOGIN>/uploads/`
  REST endpoint with a Bearer API key, capturing the project URL the
  server returns, and optionally binding a previously saved layout so
  the file opens with a pre-built arrangement of viewers and filters.
  Output is a curl/bash one-liner (or equivalent HTTP call), never a
  Datagrok package function — in-package TypeScript should use
  `grok.dapi.files.writeAsText` instead (`DG-FACT-433`).
  Use when asked to "push a csv from a shell script to a workspace",
  "curl a dataset into Datagrok from a pipeline", or "open incoming
  files with a saved viewer arrangement from outside the platform".
triggers:
  - curl a csv into a datagrok workspace
  - post a csv from a shell script
  - shareable project url from a bash batch job
  - apply saved layout to a csv uploaded over http
  - drop a file under a user account from an external pipeline
allowed-tools:
  - Read
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# upload-data

## When to use

An external process (cron job, CI pipeline, lab instrument adapter,
standalone Python notebook) has just produced a CSV and you want it
to land in Datagrok as a shareable project — ideally already wrapped
in a known viewer layout. The pipeline calls one HTTP endpoint and
gets a URL it can drop into Slack or an email.

## When NOT to use

For in-package TypeScript, use `grok.dapi.files.writeAsText(path,
csv)` / `grok.dapi.files.write(path, bytes)` — do NOT wrap the REST
endpoint in `fetch()` (`DG-FACT-433`).

## Prerequisites

- A reachable Datagrok server (`$DATAGROK_URI`, e.g.
  `https://public.datagrok.ai` or `http://localhost:8080`).
- An **API key** copied from the user profile → "API Key" link
  (`$DATAGROK_URI/u`). The REST API expects it prefixed with `Bearer `
  in the `Authorization` header (`DG-FACT-229`).
- The target user's **login** (the one that owns the resulting
  project). The article's example uses `alex.aprm`
  (`help/develop/how-to/data/upload-data.md:34`).
- A saved layout (only for the `?layout=` step) — see the related
  `layouts` skill for the JS path; the article also points to the
  manual `View | Layout | Save to Gallery` flow.

## Steps

1. **Construct the upload URL.** Fixed grammar:
   `/api/files/<LOGIN>/uploads/<your/free/path>/<file>.csv` — `<LOGIN>`
   owns the resulting project, the literal `uploads/` segment is
   required (`DG-FACT-431`).
   ```bash
   DATAGROK_URI="http://localhost:8080"
   LOGIN="alex.aprm"
   UPLOAD_URL="$DATAGROK_URI/api/files/$LOGIN/uploads/from_excel/test.csv"
   ```

2. **POST the CSV body and capture the returned URL.**
   Use `--data-binary @<path>` to preserve newlines. The response body
   is the project URL.
   ```bash
   RESPONSE=$(curl -sS --location --request POST "$UPLOAD_URL" \
     --header "Authorization: Bearer $DATAGROK_API_KEY" \
     --header "Content-Type: text/csv" \
     --data-binary @./pipeline-output.csv \
     --write-out $'\n%{http_code}')
   HTTP_CODE="${RESPONSE##*$'\n'}"
   PROJECT_URL="${RESPONSE%$'\n'*}"
   [ "$HTTP_CODE" = "200" ] || { echo "upload failed: HTTP $HTTP_CODE" >&2; exit 1; }
   case "$PROJECT_URL" in http*) ;; *) echo "unexpected body: $PROJECT_URL" >&2; exit 1;; esac
   echo "Open in Datagrok: $PROJECT_URL"
   ```
   Re-running silently overwrites the same target (`DG-FACT-432`).

3. **Bind a saved layout (optional).** Append
   `?layout=<id-or-fully-qualified-name>`. Find the ID/name in
   **Manage → Layouts → Context Panel → Links...**
   ```bash
   LAYOUT="alexaprm.superlayout"   # owner.layoutName, or the layout UUID
   curl --location --request POST "$UPLOAD_URL?layout=$LAYOUT" \
     --header "Authorization: Bearer $DATAGROK_API_KEY" \
     --header "Content-Type: text/csv" \
     --data-binary @./pipeline-output.csv
   ```

## Common failure modes

- **HTTP 401 / 403.** Add `Authorization: Bearer $DATAGROK_API_KEY`
  (the `Bearer ` prefix is mandatory — `DG-FACT-229`).
- **`?layout=mylayout` ignored.** Short name doesn't resolve across
  owners. Use the fully qualified `<owner>.<layoutname>` or the UUID.
- **CSV corruption.** Use `--data-binary @<path>`, not `--data-raw`.
- **Wrong `<LOGIN>` silently creates an unreachable upload.** Confirm
  the login from the profile page (e.g. `alex.aprm`).
- **Two pipelines collide.** Endpoint overwrites without warning
  (`DG-FACT-432`). Namespace with timestamp / run id.
- **Skill applied inside a Datagrok package.** Use
  `grok.dapi.files.writeAsText(path, df.toCsv())` instead
  (`DG-FACT-433`).

## See also

- Source: `help/develop/how-to/data/upload-data.md`,
  `help/develop/packages/rest-api.md`.
- Knowledge: `DG-FACT-229`, `431`, `432`, `433`.
- Related skills: `layouts`, `manage-credentials`, `access-data`.
