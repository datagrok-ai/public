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

This skill is for **out-of-process** callers only. If you are writing
TypeScript inside a Datagrok package (`src/package.ts`, a viewer, an
app), **do not** wrap the REST endpoint in a `fetch()` call — that
re-implements the platform's own file API, requires the user to paste
an API key, and ties the package to a hard-coded host URL. Use
`await grok.dapi.files.writeAsText(path, csv)` for text or
`await grok.dapi.files.write(path, bytes)` for binary
(`number[]`) — those are the canonical in-package paths
(`js-api/src/dapi.ts:1457,1465`), exercised by `ApiTests/src/dapi/files.ts:14,27,40`
and used by HitTriage, FileEditors, and others (`DG-FACT-433`). See
the `access-data` skill for the full `grok.dapi.files.*` surface.

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

1. **Construct the upload URL.**
   The path grammar is fixed:
   `/api/files/<LOGIN>/uploads/<your/free/path>/<file>.csv` —
   `<LOGIN>` owns the resulting project, the literal `uploads/`
   segment is required, and anything between `/uploads/` and the
   filename is a free-form path you choose (`DG-FACT-431`;
   `help/develop/how-to/data/upload-data.md:5-15`).
   ```bash
   DATAGROK_URI="http://localhost:8080"
   LOGIN="alex.aprm"
   UPLOAD_URL="$DATAGROK_URI/api/files/$LOGIN/uploads/from_excel/test.csv"
   ```
   Expected: `$UPLOAD_URL` matches the shape
   `…/api/files/alex.aprm/uploads/from_excel/test.csv`.

2. **POST the CSV body and capture the returned URL.**
   Send the file bytes with `Content-Type: text/csv`. The article's
   canonical example uses `curl --data-raw` with an inline CSV; for a
   pipeline, use `--data-binary @<path>` so newlines and special
   characters survive intact
   (`help/develop/how-to/data/upload-data.md:33-38`). One POST does
   both jobs — the response body is the URL of the uploaded project,
   so capture stdout into `$PROJECT_URL`; ask curl for the HTTP status
   on a separate line so you can fail fast without re-uploading
   (`help/develop/how-to/data/upload-data.md:17`).
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
   Expected: `$HTTP_CODE` is `200`; `$PROJECT_URL` starts with `http`
   and, when opened, renders the freshly-uploaded TableView for any
   user with read access to `$LOGIN`'s uploads. Re-running this step
   silently overwrites the same target — see failure modes below.

3. **Bind a saved layout (optional).**
   Append `?layout=<id-or-fully-qualified-name>` to apply a saved
   layout to the dataset. Find the ID or qualified name in
   **Manage → Layouts**, select the layout, open the **Context
   Panel**, hit **Links...**, copy either field; rename if ambiguous
   (`help/develop/how-to/data/upload-data.md:19-29`):
   ```bash
   LAYOUT="alexaprm.superlayout"   # owner.layoutName, or the layout UUID
   curl --location --request POST "$UPLOAD_URL?layout=$LAYOUT" \
     --header "Authorization: Bearer $DATAGROK_API_KEY" \
     --header "Content-Type: text/csv" \
     --data-binary @./pipeline-output.csv
   ```
   Expected: opening the returned URL renders the CSV with the
   layout's viewer arrangement and filter state applied. For the
   matching rules that decide which layouts the platform considers
   applicable to a given frame, see the `layouts` skill
   (`getApplicable`, `layout-id` tag).

## Common failure modes

- **HTTP 401 / 403 from the upload.** Article's curl example omits
  auth because it's a localhost demo, but any real server rejects
  unauthenticated writes. Fix: add
  `-H "Authorization: Bearer $DATAGROK_API_KEY"` (the `Bearer ` prefix
  is mandatory — `DG-FACT-229`).
- **`?layout=mylayout` is ignored / a generic table view opens.**
  Short layout name without the owner prefix doesn't resolve when
  another user owns it. Fix: copy the **fully qualified name**
  (`<owner>.<layoutname>`) or the **UUID** from **Manage → Layouts →
  Context Panel → Links...**
  (`help/develop/how-to/data/upload-data.md:24-27`).
- **CSV uploads but newlines / quotes look corrupted in the viewer.**
  `--data-raw` and shells eat backslashes and squash whitespace. Fix:
  use `--data-binary @<path>` with the file on disk (step 2 example).
- **Wrong `<LOGIN>` silently creates an unreachable upload.** The
  endpoint accepts any string for `<LOGIN>`; a typo lands the file
  under a path nobody can browse to. Fix: confirm the login exactly
  as it appears on the user's profile page, including dots (e.g.
  `alex.aprm` — `help/develop/how-to/data/upload-data.md:34`).
- **Two pipelines collide on the same `/uploads/<path>/<file>`
  target.** The REST API explicitly overwrites without warning
  (`DG-FACT-432`). Fix: namespace the path with a timestamp or run
  id (`/uploads/nightly/$(date +%Y-%m-%d)/assay.csv`).
- **Skill applied inside a Datagrok package (`fetch()` wrapper added
  to `src/package.ts`).** The REST upload surface is for out-of-process
  callers; wrapping it in TS re-implements `grok.dapi.files.*` and
  needlessly demands an API key from session code that already has
  one. Fix: drop the package function and call
  `await grok.dapi.files.writeAsText(path, df.toCsv())` for text or
  `await grok.dapi.files.write(path, bytes)` for binary
  (`js-api/src/dapi.ts:1457,1465`) — see the `access-data` skill
  (`DG-FACT-433`).

## See also

- Source articles:
  - `help/develop/how-to/data/upload-data.md`
  - `help/develop/packages/rest-api.md` (auth header, file/table
    endpoints, "replace content on re-upload")
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` —
    `DG-FACT-431` (legacy `/api/files/<LOGIN>/uploads/` path grammar),
    `DG-FACT-432` (silent overwrite on re-POST),
    `DG-FACT-229` (REST auth: `Bearer ` prefix),
    `DG-FACT-433` (in-package callers use `grok.dapi.files.writeAsText`).
- Related skills:
  - `layouts` — save / fetch `.layout` files and the
    `getApplicable` / `layout-id` tag flow referenced by step 3.
  - `manage-credentials` — pattern for storing the API key out of
    source.
  - `access-data` — fetching the uploaded CSV back from a package via
    `grok.dapi.files.*`.
