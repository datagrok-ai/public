#!/usr/bin/env bash
# Chem's engine settings on a stand, to run whole suites in Crux mode (Chem's package tests, its BDD features, the
# tests of the packages that call Chem): the package settings live on the server, so every page that starts after a
# `set` runs Chem with them.
#
#   engine-setting.sh show                      Chem's settings, per group that has some
#   engine-setting.sh set <RDKit|Crux> [<RDKit|Crux>]
#                                               ChemicalEngine (and SubstructureSearchEngine, when given) for All
#                                               users and every group that has Chem settings; what was there is saved
#   engine-setting.sh restore [<settings.json>] puts back what the last `set` saved, or the per-group maps of the file
#                                               (the properties route's JSON: {"<group id>": {...}, ...})
#
# The server keeps a group's settings as one map, replaced by each save, and shows a setting the map lacks with the
# published version's default (packages_service.dart, getPropertyValues): so `restore` stores what `show` showed, and
# an engine setting a group had unset comes back stored with its default value.
#
# HOST     the grok CLI's server alias (default localhost); the CLI must be configured for it (grok config)
# STATE    where `set` saves the settings it replaced (default ~/.grok/chem-engine-setting.json)
#
# For example, Chem's package tests in Crux mode:
#   bash engine-setting.sh set Crux && grok test --host localhost --skip-build --skip-publish; bash engine-setting.sh restore
set -euo pipefail
HOST=${HOST:-localhost}
STATE=${STATE:-$HOME/.grok/chem-engine-setting.json}
ALL_USERS=a4b45840-9a50-11e6-9cc9-8546b8bf62e6
cmd=${1:-show}

# reads JSON from stdin (a string holding JSON too: the properties route answers text) and prints the expression of it, v
json() { node -e "let s='';process.stdin.on('data',(d)=>s+=d).on('end',()=>{let v=JSON.parse(s);if(typeof v==='string')v=JSON.parse(v);console.log($1)})"; }

id=$(grok s packages list --host "$HOST" --filter 'name = "Chem"' --output json | json 'v[0].id')
current=$(grok s raw GET "/packages/$id/properties?result=false" --host "$HOST" --output json | json 'JSON.stringify(v)')

post() {  # post <group id> <settings JSON>
  grok s raw POST "/packages/$id/properties?groupId=$1" --data "$2" --host "$HOST" --output json > /dev/null
}

case "$cmd" in
show)
  echo "$current" | json 'Object.entries(v).map(([g, p]) => `${g}${g === "'$ALL_USERS'" ? " (All users)" : ""}: ` +
    `ChemicalEngine ${p.ChemicalEngine ?? "(default RDKit)"}, SubstructureSearchEngine ${p.SubstructureSearchEngine ?? "(default)"}`).join("\n")'
  ;;
set)
  engine=${2:?usage: engine-setting.sh set <RDKit|Crux> [<RDKit|Crux>]}
  search=${3:-}
  for e in "$engine" ${search:+"$search"}; do
    case "$e" in
      RDKit|Crux) ;;
      *) echo "engine-setting.sh: '$e' is not an engine (RDKit or Crux; the names are case-sensitive)" >&2; exit 2 ;;
    esac
  done
  mkdir -p "$(dirname "$STATE")"
  [ -f "$STATE" ] || echo "$current" > "$STATE"   # a second `set` keeps the first one's original values
  for g in $(echo "$current" | json "[...new Set([...Object.keys(v), '$ALL_USERS'])].join(' ')"); do
    settings=$(echo "$current" | json "JSON.stringify({...(v['$g'] ?? {}), ChemicalEngine: '$engine'$( [ -n "$search" ] && echo ", SubstructureSearchEngine: '$search'")})")
    post "$g" "$settings"
  done
  echo "Chem: ChemicalEngine $engine${search:+, SubstructureSearchEngine $search}; the replaced settings are in $STATE"
  ;;
restore)
  from=${2:-$STATE}
  [ -f "$from" ] || { echo "nothing to restore ($from does not exist)"; exit 0; }
  saved=$(json 'JSON.stringify(v)' < "$from")
  # the groups `set` wrote: those it saved, and All users (none saved: it had none)
  for g in $(echo "$saved" | json "[...new Set([...Object.keys(v), '$ALL_USERS'])].join(' ')"); do
    post "$g" "$(echo "$saved" | json "JSON.stringify(v['$g'] ?? {})")"
  done
  # a server that merges the posted settings into the stored ones keeps a key the saved ones lack: the engine settings
  # go back to their defaults then
  after=$(grok s raw GET "/packages/$id/properties?result=false" --host "$HOST" --output json | json 'JSON.stringify(v)')
  for g in $(echo "$after" | json "Object.keys(v).join(' ')"); do
    fix=$(node -e "const a = JSON.parse(process.argv[1])['$g'] ?? {}, o = JSON.parse(process.argv[2])['$g'] ?? {};
      const out = {...a}; let changed = false;
      if (!('ChemicalEngine' in o) && a.ChemicalEngine && a.ChemicalEngine !== 'RDKit') { out.ChemicalEngine = 'RDKit'; changed = true; }
      if (!('SubstructureSearchEngine' in o) && a.SubstructureSearchEngine && a.SubstructureSearchEngine !== 'Crux') { out.SubstructureSearchEngine = 'Crux'; changed = true; }
      console.log(changed ? JSON.stringify(out) : '')" "$after" "$saved")
    if [ -n "$fix" ]; then
      post "$g" "$fix"
      echo "the server kept a setting $g did not have: set to its default"
    fi
  done
  [ "$from" = "$STATE" ] && rm -f "$STATE"
  echo "Chem: engine settings restored from $from"
  ;;
*)
  echo "usage: engine-setting.sh show | set <RDKit|Crux> [<RDKit|Crux>] | restore [<settings.json>]" >&2
  exit 2
  ;;
esac
