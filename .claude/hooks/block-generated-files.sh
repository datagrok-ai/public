#!/bin/bash
# PreToolUse hook: block editing auto-generated .g.ts and .api.g.ts files
input=$(cat)
file=$(echo "$input" | sed -n 's/.*"file_path"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' | head -1)

if echo "$file" | grep -qE '\.(g\.ts|api\.g\.ts)$'; then
  echo '{"hookSpecificOutput": {"hookEventName": "PreToolUse", "permissionDecision": "deny", "permissionDecisionReason": "This is an auto-generated file (.g.ts / .api.g.ts). Run grok api to regenerate instead of editing manually."}}'
fi
