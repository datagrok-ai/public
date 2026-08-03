#!/usr/bin/env bash
set -euo pipefail

# npm dependency audit of ONE just-published Datagrok npm package → OpenVEX +
# CSV + HTML + meta on the public VEX surface (GROK-18695). Mirrors the weekly
# Jenkins Vulnerability-Scan npm pass (infra/jenkins/_scripts/vuln_scan/
# npm-audit-vex.sh in the reddata repo) for use in GitHub Actions at package
# publish time, so the package's npm VEX is available the moment it is
# published (the weekly job then skips the version via the S3 delta).
#
# Audits the PRODUCTION dependency tree only (--omit=dev): devDependencies
# never reach a consumer install. The workspace package-lock.json is used, so
# the audit reflects the tree the shipped webpack bundle was built from
# (including the package's own `overrides`).
#
# Requires on PATH: npm, jq, aws. Env: AWS_ACCESS_KEY_ID/SECRET (S3 upload).
# Usage: npm-audit-vex.sh <package_dir>   e.g. npm-audit-vex.sh packages/Chem

PKG_DIR="${1:?package dir required}"
[[ -f "${PKG_DIR}/package.json" ]] || { echo "No package.json in ${PKG_DIR}"; exit 1; }

S3_BUCKET="${S3_BUCKET:-datagrok-data}"
S3_PREFIX="${S3_PREFIX:-vex}"
AWS_REGION="${AWS_REGION:-us-east-2}"
PUBLIC_BASE="${PUBLIC_BASE:-https://data.datagrok.ai/${S3_PREFIX}}"
RUN_TS="$(date -u +%FT%TZ)"

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing $1"; exit 1; }; }
need npm; need jq; need aws

NAME="$(jq -r '.name' "${PKG_DIR}/package.json")"
VER="$(jq -r '.version' "${PKG_DIR}/package.json")"
SAFE="${NAME#@}"; SAFE="${SAFE//\//-}"      # @datagrok/chem → datagrok-chem
PURL="pkg:npm/${NAME/@/%40}@${VER}"

cd "$PKG_DIR"
[[ -f package-lock.json ]] || \
  npm install --package-lock-only --ignore-scripts --no-audit --no-fund --loglevel=error

audit="$(mktemp)"
npm audit --omit=dev --package-lock-only --json > "$audit" 2>/dev/null || true
jq -e '.auditReportVersion // .metadata' "$audit" >/dev/null 2>&1 \
  || { echo "npm audit produced no report for ${NAME}@${VER}"; exit 1; }

# Flatten to one finding per (advisory, dependency); advisory objects live in
# the .via arrays of the packages they are about, string entries are
# transitive references. Installed versions come from the lockfile.
findings="$(mktemp)"
jq --slurpfile lock package-lock.json '
  def sevmap: {"critical":"CRITICAL","high":"HIGH","moderate":"MEDIUM",
               "low":"LOW","info":"MINIMAL"}[.] // "UNKNOWN";
  [ (.vulnerabilities // {}) | to_entries[] | .value as $v | ($v.via // [])[]
    | select(type == "object")
    | { id: (if ((.url // "") | test("/advisories/|/vulnerability/"))
               then ((.url) | split("/") | last)
             elif (.source != null) then "NPM-ADVISORY-\(.source)"
             else empty end),
        pkg:   (.name // $v.name),
        sev:   ((.severity // "unknown") | sevmap),
        cvss:  ((.cvss.score // "") | tostring | if . == "0" then "" else . end),
        range: (.range // ""),
        desc:  (.title // ""),
        url:   (.url // ""),
        fixed: (if ($v.fixAvailable | type) == "object"
                  then "\($v.fixAvailable.name)@\($v.fixAvailable.version)"
                elif $v.fixAvailable == true then "available" else "" end),
        installed: ($lock[0].packages[($v.nodes // [])[0] // ""].version // "") } ]
  | unique_by([.id, .pkg])' "$audit" > "$findings"

read -r C H M L T <<< "$(jq -r '
  unique_by(.id) as $f
  | [ ($f | map(select(.sev=="CRITICAL")) | length),
      ($f | map(select(.sev=="HIGH"))     | length),
      ($f | map(select(.sev=="MEDIUM"))   | length),
      ($f | map(select(.sev=="LOW" or .sev=="MINIMAL")) | length),
      ($f | length) ] | @tsv' "$findings" | tr '\t' ' ')"
echo "→ ${NAME}@${VER}: ${C} critical / ${H} high / ${M} medium / ${L} low (${T} advisories)"

vex="$(mktemp)"; csv="$(mktemp)"; html="$(mktemp)"; meta="$(mktemp)"

# OpenVEX 0.2.0 — one statement per distinct advisory.
jq --arg id "${PUBLIC_BASE}/npm/${SAFE}/${VER}.json" --arg ts "$RUN_TS" \
   --arg product "$PURL" '
  group_by(.id)
  | map({ vulnerability: { name: .[0].id },
          products: [ { "@id": $product } ],
          status: "affected" }
        + ( if (map(.fixed) | any(. != "")) then
              { action_statement: "A fix is available in a newer dependency version; upgrade the affected npm dependency." }
            else {} end ))
  | { "@context": "https://openvex.dev/ns/v0.2.0",
      "@id": $id, author: "Datagrok", timestamp: $ts, version: 1,
      statements: . }' "$findings" > "$vex"

# CSV — same column layout as the docker image CSVs.
{
  echo 'repo,tag,package,cve,severity,cvss,installed_version,fixed_version,status'
  jq -r --arg repo "$NAME" --arg tag "$VER" '
    .[] | [ $repo, $tag, .pkg, .id, .sev, .cvss, .installed, .fixed, "affected" ] | @csv' "$findings"
} > "$csv"

# Human-readable HTML report (same chrome as the rest of the VEX surface).
DOCS_URL="${DOCS_URL:-https://datagrok.ai/help/datagrok/solutions/teams/it/security}"
badge() { if [[ "${1:-0}" -gt 0 ]]; then echo "<span class=\"badge $2\">$1</span>"; else echo "<span class=\"badge b-zero\">0</span>"; fi; }
{
  cat <<HTML
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>${NAME}@${VER} — npm dependency audit</title>
<style>
  :root { color-scheme: light dark; }
  body { font: 14px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
         margin: 0; padding: 24px 32px; color: #1a1a1a; background: #fff; }
  h1 { font-size: 22px; margin: 0 0 4px; }
  .sub { color: #666; margin: 0 0 16px; }
  nav { margin: 0 0 20px; font-size: 13px; }
  nav a { margin-right: 14px; }
  a { color: #1f6feb; text-decoration: none; }
  a:hover { text-decoration: underline; }
  table { border-collapse: collapse; width: 100%; margin: 8px 0 4px; }
  th, td { text-align: left; padding: 6px 10px; border-bottom: 1px solid #e2e2e2;
           vertical-align: top; }
  th { background: #f6f8fa; font-weight: 600; white-space: nowrap; }
  td.num { text-align: right; font-variant-numeric: tabular-nums; }
  .badge { display: inline-block; min-width: 20px; text-align: center; padding: 1px 7px;
           border-radius: 10px; font-size: 12px; font-weight: 600; color: #fff; }
  .b-crit { background: #b31d28; } .b-high { background: #d1541f; }
  .b-med  { background: #b08800; } .b-low  { background: #57606a; }
  .b-zero { background: #dfe3e8; color: #57606a; }
  tr.sev-critical td:first-child { border-left: 3px solid #b31d28; }
  tr.sev-high     td:first-child { border-left: 3px solid #d1541f; }
  tr.sev-medium   td:first-child { border-left: 3px solid #b08800; }
  tr.sev-low, tr.sev-minimal td:first-child { border-left: 3px solid #57606a; }
  footer { margin: 32px 0 8px; color: #888; font-size: 12px; }
  code { background: #f0f1f3; padding: 1px 4px; border-radius: 4px; }
</style>
</head>
<body>
<h1>${NAME}@${VER} — npm dependency audit</h1>
<p class="sub">$(badge "$C" b-crit) critical · $(badge "$H" b-high) high · $(badge "$M" b-med) medium · $(badge "$L" b-low) low · <b>${T}</b> total (deduped by advisory)</p>
<nav>
  <a href="../../index.html">&larr; Vulnerability index</a>
  <a href="./${VER}.json">OpenVEX</a>
  <a href="./${VER}.csv">CSV</a>
  <a href="../../npm.json">Composite (npm)</a>
  <a href="${DOCS_URL}">Docs</a>
</nav>
<p class="sub">NPM package <code>${NAME}@${VER}</code> — production dependency tree audited with <code>npm audit</code>.</p>
<table>
<thead><tr>
  <th>Advisory</th><th>Severity</th><th>Dependency</th>
  <th>Installed</th><th>Vulnerable range</th><th>Fixed by</th><th>CVSS</th><th>Description</th>
</tr></thead>
<tbody>
HTML
  jq -r '
    def esc: tostring
      | gsub("&";"&amp;") | gsub("<";"&lt;") | gsub(">";"&gt;") | gsub("\"";"&quot;");
    def sevrank: {"CRITICAL":0,"HIGH":1,"MEDIUM":2,"LOW":3,"MINIMAL":4,"UNKNOWN":5}[.] // 6;
    sort_by([(.sev|sevrank), .id]) | .[]
    | "<tr class=\"sev-\(.sev|ascii_downcase)\">"
      + "<td><a href=\"\(.url|esc)\">\(.id|esc)</a></td>"
      + "<td>\(.sev|esc)</td>"
      + "<td>\(.pkg|esc)</td>"
      + "<td>\(if .installed=="" then "&mdash;" else (.installed|esc) end)</td>"
      + "<td>\(.range|esc)</td>"
      + "<td>\(if .fixed=="" then "&mdash;" else (.fixed|esc) end)</td>"
      + "<td class=\"num\">\(.cvss|esc)</td>"
      + "<td>\(.desc|esc)</td></tr>"' "$findings"
  cat <<HTML
</tbody>
</table>
HTML
  [[ "$T" -eq 0 ]] && echo '<p class="sub">No known vulnerabilities in the production dependency tree.</p>'
  cat <<HTML
<footer>
  Generated ${RUN_TS} at package publish time (GitHub Actions, GROK-18695).
  Scanner: npm audit (production dependencies). Format: <a href="https://openvex.dev/">OpenVEX 0.2.0</a>.
  See the <a href="${DOCS_URL}">security documentation</a> for the remediation policy.
</footer>
</body>
</html>
HTML
} > "$html"

# Meta record feeding the index (rebuilt by the Jenkins REINDEX_ONLY trigger).
jq -n --arg repo "npm/${SAFE}" --arg name "$NAME" --arg tag "$VER" --arg ts "$RUN_TS" \
   --argjson c "$C" --argjson h "$H" --argjson m "$M" --argjson l "$L" --argjson t "$T" '
  { repo: $repo, name: $name, tag: $tag, digest: "", category: "npm",
    description: "Production npm dependency audit",
    critical: $c, high: $h, medium: $m, low: $l, total: $t, scanned: $ts }' > "$meta"

put() { aws s3 cp "$1" "s3://${S3_BUCKET}/${S3_PREFIX}/$2" --acl public-read --content-type "$3" --region "$AWS_REGION" >/dev/null; }
put_mutable() { aws s3 cp "$1" "s3://${S3_BUCKET}/${S3_PREFIX}/$2" --acl public-read --content-type "$3" --region "$AWS_REGION" --cache-control "public, max-age=300" >/dev/null; }
# At publish time this version IS the newest, so refresh latest + meta too.
put "$vex"  "npm/${SAFE}/${VER}.json" application/json
put "$csv"  "npm/${SAFE}/${VER}.csv"  text/csv
put "$html" "npm/${SAFE}/${VER}.html" text/html
put_mutable "$vex"  "npm/${SAFE}/latest.json"      application/json
put_mutable "$csv"  "npm/${SAFE}/latest.csv"       text/csv
put_mutable "$html" "npm/${SAFE}/latest.html"      text/html
put_mutable "$meta" "npm/${SAFE}/latest.meta.json" application/json
echo "✔ Published ${PUBLIC_BASE}/npm/${SAFE}/${VER}.json (+ .csv, .html, latest.*)"
