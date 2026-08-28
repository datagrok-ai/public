#!/usr/bin/env bash
# Enroll Datagrok npm packages in npm trusted publishing (OIDC).
#
# Maps every publishable package in this repo to the workflow that publishes it
# and registers it with `npm trust github`. Idempotent: packages that already
# have a trusted publisher, are not on the registry yet, or are not ours are
# skipped.
#
# Usage:
#   ./.github/scripts/npm-trust-enroll.sh                 # everything
#   ./.github/scripts/npm-trust-enroll.sh --dry-run       # show what would happen
#   ./.github/scripts/npm-trust-enroll.sh @datagrok/chem  # one package
#
# Requires npm >= 11.15.0, `npm login` as a user with write access to the
# @datagrok / @datagrok-libraries / @datagrok-misc orgs, and account-level 2FA.
# See .github/NPM_TRUSTED_PUBLISHING.md.
set -uo pipefail

REPO='datagrok-ai/public'
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

DRY_RUN=''
ONLY=''
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN='--dry-run' ;;
    -*) echo "unknown option: $arg" >&2; exit 2 ;;
    *) ONLY="$arg" ;;
  esac
done

npm_major_minor="$(npm --version | awk -F. '{printf "%d%03d", $1, $2}')"
if [ "$npm_major_minor" -lt 11015 ]; then
  echo "npm $(npm --version) is too old; 'npm trust' needs >= 11.15.0. Run: npm install -g npm@latest" >&2
  exit 1
fi

enrolled=0 skipped=0 failed=0

enroll() {
  local dir="$1" workflow="$2"
  [ -f "$dir/package.json" ] || return 0

  local name
  name="$(jq -r '.name // empty' "$dir/package.json")"
  [ -n "$name" ] || return 0
  [ -z "$ONLY" ] || [ "$ONLY" = "$name" ] || return 0

  case "$name" in
    @datagrok/*|@datagrok-libraries/*|@datagrok-misc/*|datagrok-api|datagrok-tools) ;;
    *) echo "SKIP  $name — not a Datagrok-owned name ($dir)"; skipped=$((skipped + 1)); return 0 ;;
  esac

  if [ "$(jq -r '.private // false' "$dir/package.json")" = 'true' ]; then
    echo "SKIP  $name — private"; skipped=$((skipped + 1)); return 0
  fi

  local escaped code
  escaped="${name//\//%2f}"
  code="$(curl -s -o /dev/null -w '%{http_code}' --retry 2 "https://registry.npmjs.org/${escaped}")"
  if [ "$code" != '200' ]; then
    echo "SKIP  $name — not on the registry yet (publish its first version by hand, then re-run)"
    skipped=$((skipped + 1)); return 0
  fi

  # One call per package: npm trust is OTP-gated, and the "skip 2FA" window is
  # only 5 minutes. Registering an already-registered package is a plain error,
  # so classify the failure rather than paying for a separate `npm trust list`.
  local out rc err_code err_text
  # shellcheck disable=SC2086
  # stderr is deliberately not redirected: npm prints the OTP approval URL there
  # and waits, so swallowing it would look like a silent hang.
  out="$(npm trust github "$name" --repo "$REPO" --file "$workflow" \
    --allow-publish --yes --json $DRY_RUN)"
  rc=$?
  err_code="$(jq -r '.error.code // empty' 2>/dev/null <<<"$out")"
  err_text="$(jq -r '.error.summary // .error.detail // empty' 2>/dev/null <<<"$out")"
  auth_url="$(jq -r '.error.authUrl // empty' 2>/dev/null <<<"$out")"
  # npm masks the UUID in its human-readable output but not in --json, so
  # recover the approval URL even when the payload will not parse.
  [ -n "$auth_url" ] || auth_url="$(grep -oE 'https://www\.npmjs\.com/auth/cli/[0-9a-fA-F-]{36}' <<<"$out" | head -1)"
  [ -n "$err_code" ] || err_code="$(sed -n 's/.*"code"[[:space:]]*:[[:space:]]*"\([A-Z]*\)".*/\1/p' <<<"$out" | head -1)"

  if [ "$rc" -eq 0 ]; then
    echo "TRUST $name -> $workflow"
    enrolled=$((enrolled + 1))
  elif [ "$err_code" = 'EOTP' ]; then
    echo >&2
    echo "npm requires two-factor approval before it will accept trust changes." >&2
    echo "Open this URL, tick 'skip two-factor authentication for the next" >&2
    echo "5 minutes', approve, then re-run this script:" >&2
    echo >&2
    echo "    $auth_url" >&2
    echo >&2
    echo "(npm hides this UUID in its own error output; it is read from --json.)" >&2
    echo "The URL is single-use and expires in a few minutes." >&2
    exit 3
  elif [ -n "$err_text" ] && grep -qiE 'already|exists|conflict' <<<"$err_text"; then
    echo "SKIP  $name — already has a trusted publisher"; skipped=$((skipped + 1))
  else
    echo "FAIL  $name — ${err_code:-exit $rc}: ${err_text:-see npm output above}" >&2
    failed=$((failed + 1))
    # A failure before anything has succeeded is systemic (auth, network,
    # permissions) — bail rather than replay it across every package.
    if [ "$enrolled" -eq 0 ]; then
      echo >&2
      echo "Stopping: first package failed and nothing has succeeded yet." >&2
      echo "If npm printed an auth URL above, open it, tick 'skip two-factor" >&2
      echo "authentication for the next 5 minutes', then re-run this script." >&2
      exit 3
    fi
  fi
  # Stay under the registry rate limit; pointless when nothing is being written.
  [ -n "$DRY_RUN" ] || sleep 2
}

enroll "$ROOT/js-api" 'js-api.yml'
enroll "$ROOT/tools" 'tools.yml'
for d in "$ROOT"/libraries/*/; do enroll "$d" 'libraries.yaml'; done
for d in "$ROOT"/misc/*/; do enroll "$d" 'misc.yaml'; done
for d in "$ROOT"/packages/*/; do enroll "$d" 'packages.yaml'; done

echo
echo "enrolled=$enrolled skipped=$skipped failed=$failed"

if [ -n "$ONLY" ] && [ $((enrolled + skipped + failed)) -eq 0 ]; then
  echo "No package matched '$ONLY'." >&2
  exit 2
fi
[ "$failed" -eq 0 ]
