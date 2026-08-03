#!/bin/bash
# Tests for the install scripts' release resolution.
#
# A compose file and the images it starts must come from the same release. The scripts work
# that out by reading the `branch` label off the image they are about to run, so these cases
# stub `docker` and assert the ref/release they derive, and that a mismatched compose file on
# disk is refused. No images are pulled and nothing is started.
#
#   bash test_install_scripts.sh

set -o nounset

HERE="$(cd -- "$(dirname -- "$0")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
pass=0; fail=0

# `docker` stub: reports the label named by STUB_BRANCH and swallows everything else.
mkdir -p "$WORK/bin"
cat > "$WORK/bin/docker" <<'EOF'
#!/bin/bash
# STUB_BRANCH is what the registry reports, STUB_LOCAL_BRANCH what the copy on disk reports.
# Unset STUB_LOCAL_BRANCH defaults to the registry answer (an up-to-date install).
case "${1:-} ${2:-}" in
  "buildx imagetools")
    [ -n "${STUB_BRANCH:-}" ] && echo "$STUB_BRANCH"
    exit 0 ;;
  "image inspect")
    echo "${STUB_LOCAL_BRANCH-${STUB_BRANCH:-}}"
    exit 0 ;;
  "images -q") echo "stubimageid"; exit 0 ;;
  "info ")     exit 0 ;;
esac
[ "${1:-}" = "info" ] && exit 0
exit 0
EOF
chmod +x "$WORK/bin/docker"
for b in xdg-open open; do printf '#!/bin/bash\nexit 0\n' > "$WORK/bin/$b"; chmod +x "$WORK/bin/$b"; done
export PATH="$WORK/bin:$PATH"

ok() { printf 'ok   %s\n' "$1"; pass=$((pass+1)); }
no() { printf 'FAIL %s\n' "$1"; fail=$((fail+1)); }

setup() { # setup <script>
  rm -rf "$WORK/run"; mkdir -p "$WORK/run"
  cp "$HERE/$1" "$WORK/run/"
  sed -i 's/^timeout=30/timeout=1/' "$WORK/run/$1"
}

resolves() { # resolves <script> <label> <tag> <want-release> <want-ref>
  setup "$1"
  local out
  out=$(cd "$WORK/run" && STUB_BRANCH="$2" DATAGROK_VERSION="$3" bash "$1" version 2>&1 | tail -1)
  if [ "$out" = "$3 -> $4 ($5)" ]; then
    ok "$3 with label '${2:-<none>}' -> $4 ($5)"
  else
    no "$3 with label '${2:-<none>}': got '$out', wanted '$3 -> $4 ($5)'"
  fi
}

refuses() { # refuses <script> <label> <tag> <desc>
  setup "$1"
  ( cd "$WORK/run" && STUB_BRANCH="$2" DATAGROK_VERSION="$3" bash "$1" version >/dev/null 2>&1 )
  [ "$?" = "255" ] && ok "refuses $4" || no "refuses $4 (exit $?)"
}

marker_case() { # marker_case <script> <compose-marker> <label> <tag> <want-exit> <desc>
  setup "$1"
  local cf
  cf=$(grep -oE 'compose_config_name="[^"]+"' "$WORK/run/$1" | sed 's/.*="//;s/"//')
  if [ -n "$2" ]; then
    printf 'x-datagrok-release: %s\nservices: {}\n' "$2" > "$WORK/run/$cf"
  else
    printf 'services: {}\n' > "$WORK/run/$cf"
  fi
  ( cd "$WORK/run" && STUB_BRANCH="$3" DATAGROK_VERSION="$4" bash "$1" start >/dev/null 2>&1 )
  local rc=$?
  [ "$rc" = "$5" ] && ok "$6" || no "$6 (exit $rc, wanted $5)"
}

for s in datagrok-install-local.sh datagrok-install-local-macos-silicon.sh; do
  echo "== $s =="

  # The label is authoritative, whatever the tag is called.
  resolves "$s" "release/1.27.7" "latest"        "1.27.7"        "release/1.27.7"
  resolves "$s" "release/1.27.7" "1.27.7"        "1.27.7"        "release/1.27.7"
  resolves "$s" "release/1.28.0" "1.28.0-rc"     "1.28.0"        "release/1.28.0"
  resolves "$s" "master"         "bleeding-edge" "bleeding-edge" "master"
  # A feature-branch build tracks master's compose file.
  resolves "$s" "GROK-20452-wip" "bleeding-edge" "bleeding-edge" "master"

  # A moving tag must report where the registry points now, not where a copy pulled months
  # ago points. Resolving from the stale local image would install the wrong release.
  stale() { # stale <script> <registry> <local> <tag> <want-release> <want-ref>
    setup "$1"
    local out
    out=$(cd "$WORK/run" && STUB_BRANCH="$2" STUB_LOCAL_BRANCH="$3" DATAGROK_VERSION="$4" \
          bash "$1" version 2>&1 | tail -1)
    [ "$out" = "$4 -> $5 ($6)" ] \
      && ok "$4: registry $2 wins over stale local $3" \
      || no "$4: registry $2 vs local $3 gave '$out', wanted '$4 -> $5 ($6)'"
  }
  stale "$s" "release/1.27.7" "release/1.27.3" "latest" "1.27.7" "release/1.27.7"
  # With the registry unreachable, an existing install still resolves from what it has.
  stale "$s" "" "release/1.27.3" "latest" "1.27.3" "release/1.27.3"

  # No label (images built before it existed): fall back to reading the tag.
  resolves "$s" "" "1.27.4"        "1.27.4"        "release/1.27.4"
  resolves "$s" "" "1.27.4-rc"     "1.27.4"        "release/1.27.4"
  resolves "$s" "" "bleeding-edge" "bleeding-edge" "master"
  refuses  "$s" "" "some-tag"      "an unlabelled image with an unreadable tag"

  # The compose file on disk must belong to the release the image resolved to.
  marker_case "$s" "1.27.7"        "release/1.27.7" "latest" 0   "matched compose starts"
  marker_case "$s" "1.27.3"        "release/1.27.7" "latest" 255 "stale compose refused"
  marker_case "$s" "bleeding-edge" "release/1.27.7" "latest" 255 "master compose vs a release refused"
  marker_case "$s" "1.27.7"        "master" "bleeding-edge"  255 "release compose vs master refused"
  marker_case "$s" ""              "release/1.27.7" "latest" 0   "unmarked legacy compose warns and runs"
done

echo
printf 'TOTAL: %d passed, %d failed\n' "$pass" "$fail"
[ "$fail" = "0" ]
