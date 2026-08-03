#!/usr/bin/env bash
set -euo pipefail

# GCP Container Analysis scan of ONE image -> OpenVEX + CSV -> public S3.
# Mirrors the Jenkins Vulnerability-Scan job (GROK-18695) for use in GitHub
# Actions at package-publish time, so a package image's VEX is available the
# moment it is published (the weekly Jenkins job then skips it via the S3 delta).
#
# Requires on PATH: gcloud, jq, aws, docker (all preinstalled on GH ubuntu runners).
# Env: GOOGLE_APPLICATION_CREDENTIALS (SA key file), AWS_ACCESS_KEY_ID/SECRET (env).
# Usage: gcp-vuln-scan.sh <image_ref> <repo> <tag>
#   e.g. gcp-vuln-scan.sh datagrok/admetica:1.3.1 admetica 1.3.1

IMAGE_REF="${1:?image_ref required}"
REPO="${2:?repo required}"
TAG="${3:?tag required}"

AR_LOCATION="${AR_LOCATION:-us-east4}"
GCP_PROJECT="${GCP_PROJECT:-axial-reference-466209-j9}"
AR_REPO="${AR_REPO:-datagrok}"
AR_HOST="${AR_LOCATION}-docker.pkg.dev"
AR_PREFIX="${AR_HOST}/${GCP_PROJECT}/${AR_REPO}"
DST_REF="${AR_PREFIX}/${REPO}:${TAG}"

S3_BUCKET="${S3_BUCKET:-datagrok-data}"
S3_PREFIX="${S3_PREFIX:-vex}"
AWS_REGION="${AWS_REGION:-us-east-2}"
PUBLIC_BASE="${PUBLIC_BASE:-https://data.datagrok.ai/${S3_PREFIX}}"
RUN_TS="$(date -u +%FT%TZ)"

SCAN_POLL_ATTEMPTS="${SCAN_POLL_ATTEMPTS:-60}"
SCAN_POLL_INTERVAL="${SCAN_POLL_INTERVAL:-15}"
SCAN_EMPTY_GRACE="${SCAN_EMPTY_GRACE:-4}"
SCAN_DESCRIBE_TIMEOUT="${SCAN_DESCRIBE_TIMEOUT:-300}"

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing $1"; exit 1; }; }
need gcloud; need jq; need aws; need docker

gcloud auth activate-service-account --key-file "${GOOGLE_APPLICATION_CREDENTIALS}" >/dev/null
gcloud config set project "$GCP_PROJECT" >/dev/null
gcloud auth configure-docker "${AR_HOST}" --quiet >/dev/null 2>&1

# Always remove the pushed copy from GAR (delete needs --delete-tags for a tagged
# image), incl. when `timeout` sends SIGTERM.
cleanup() {
  local digest
  digest="$(gcloud artifacts docker images describe "${DST_REF}" --format='get(image_summary.digest)' 2>/dev/null || true)"
  if [[ -n "$digest" ]]; then
    gcloud artifacts docker images delete "${AR_PREFIX}/${REPO}@${digest}" --delete-tags --quiet || true
  else
    gcloud artifacts docker images delete "${DST_REF}" --delete-tags --quiet || true
  fi
  docker rmi -f "${DST_REF}" >/dev/null 2>&1 || true
}
trap cleanup EXIT
trap 'exit 143' INT TERM

echo "→ Pushing ${IMAGE_REF} → ${DST_REF}"
docker image inspect "${IMAGE_REF}" >/dev/null 2>&1 || docker pull "${IMAGE_REF}"
docker tag "${IMAGE_REF}" "${DST_REF}"
docker push "${DST_REF}"

scan_json="$(mktemp)"
empty_confirms=0
ok=false
for ((i = 1; i <= SCAN_POLL_ATTEMPTS; i++)); do
  # `describe` takes the full ref (no --location flag) and its payload can be big.
  if timeout "${SCAN_DESCRIBE_TIMEOUT}" gcloud artifacts docker images describe "${DST_REF}" \
       --format=json --show-package-vulnerability > "${scan_json}.tmp" 2>/dev/null; then
    if jq -e '.discovery_summary.discovery[0].discovery.analysisStatus != "FINISHED_SUCCESS"' \
         "${scan_json}.tmp" >/dev/null 2>&1; then
      echo "  ⏳ scan in progress ($i/${SCAN_POLL_ATTEMPTS})"; rm -f "${scan_json}.tmp"; sleep "$SCAN_POLL_INTERVAL"; continue
    fi
    vc="$(jq '[.package_vulnerability_summary.vulnerabilities // {} | to_entries[] | .value[]] | length' "${scan_json}.tmp" 2>/dev/null || echo 0)"
    if [[ "$vc" -eq 0 && "$empty_confirms" -lt "$SCAN_EMPTY_GRACE" ]]; then
      empty_confirms=$((empty_confirms + 1)); echo "  ⏳ finished, awaiting vulnerability data (grace ${empty_confirms}/${SCAN_EMPTY_GRACE})"
      rm -f "${scan_json}.tmp"; sleep "$SCAN_POLL_INTERVAL"; continue
    fi
    mv "${scan_json}.tmp" "$scan_json"; ok=true; echo "  ✔ scan complete (${vc} vulnerability entries)"; break
  fi
  echo "  ⏳ waiting for report ($i/${SCAN_POLL_ATTEMPTS})"; sleep "$SCAN_POLL_INTERVAL"
done
$ok || { echo "✖ scan did not complete for ${DST_REF}; skipping"; exit 1; }

digest="$(jq -r '.image_summary.digest // ""' "$scan_json")"
vex="$(mktemp)"; csv="$(mktemp)"

# OpenVEX 0.2.0 — one statement per distinct CVE.
jq --arg id "${PUBLIC_BASE}/${REPO}/${TAG}.json" --arg ts "$RUN_TS" \
   --arg product "pkg:oci/${REPO}@${digest}?tag=${TAG}" '
  [ (.package_vulnerability_summary.vulnerabilities // {} | to_entries[] | .value[]) ]
  | map({ cve: ((.noteName // "") | split("/") | last), fix: (.vulnerability.fixAvailable // false) })
  | map(select(.cve | startswith("CVE") or startswith("GHSA") or startswith("PYSEC")))
  | group_by(.cve)
  | map({ vulnerability: { name: .[0].cve }, products: [ { "@id": $product } ], status: "affected" }
        + ( if (map(.fix) | any) then
              { action_statement: "A fix is available in a newer package version; upgrade the affected package." }
            else {} end ))
  | { "@context": "https://openvex.dev/ns/v0.2.0", "@id": $id, author: "Datagrok",
      timestamp: $ts, version: 1, statements: . }' "$scan_json" > "$vex"

# CSV — one row per affected package/CVE.
{
  echo 'repo,tag,package,cve,severity,cvss,installed_version,fixed_version,status'
  jq -r --arg repo "$REPO" --arg tag "$TAG" '
    (.package_vulnerability_summary.vulnerabilities // {} | to_entries[] | .value[])
    | . as $o | ($o.vulnerability.packageIssue // [{}])[]
    | [ $repo, $tag, (.affectedPackage // ""), (($o.noteName // "") | split("/") | last),
        ($o.vulnerability.effectiveSeverity // $o.vulnerability.severity // ""),
        (($o.vulnerability.cvssScore // "") | tostring),
        (.affectedVersion.fullName // .affectedVersion.name // ""),
        (.fixedVersion.fullName // .fixedVersion.name // ""), "affected" ] | @csv' "$scan_json"
} > "$csv"

put() { aws s3 cp "$1" "s3://${S3_BUCKET}/${S3_PREFIX}/$2" --acl public-read --content-type "$3" --region "$AWS_REGION" >/dev/null; }
# At publish time this version IS the newest, so refresh latest too.
put "$vex" "${REPO}/${TAG}.json"    application/json
put "$csv" "${REPO}/${TAG}.csv"     text/csv
put "$vex" "${REPO}/latest.json"    application/json
put "$csv" "${REPO}/latest.csv"     text/csv
echo "✔ Published ${PUBLIC_BASE}/${REPO}/${TAG}.json (+ latest, + .csv)"
