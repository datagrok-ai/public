# Grokky skill eval harness

Headless eval harness for Grokky Claude Code skills. Calls the
`@anthropic-ai/claude-agent-sdk` directly with the same prompt + skill
assembly used by `claude-runtime/src/server.ts`, then runs static rubric
checks on each response.

## Install

```bash
cd dockerfiles/claude-runtime/skills-dev/_eval-harness
npm install
```

## Run

```bash
export ANTHROPIC_API_KEY=sk-ant-...
npx tsx run.ts ../df-and-columns
```

The harness writes a timestamped report to
`<skill-dir>/eval/runs/<YYYY-MM-DD-HHMM>/report.md`.

## How it works

1. Reads `<skill-dir>/SKILL.md` and `<skill-dir>/eval/prompts.json`.
2. For each prompt, calls `query()` with:
   - `systemPrompt` = `DATAGROK_PROMPT` (verbatim copy from `server.ts`)
     + `datagrok-exec` skill + the skill under test. Other inlined skills
     are skipped to keep the eval focused on the one skill.
   - `model: 'sonnet'`, `effort: 'low'`, `thinking: disabled` (matches the
     container).
   - No tools allowed: this is a pure prompt-+-skill eval; we want to see
     what Claude **writes**, not what it does with Bash/Grep. Tool use would
     blow up costs and not test the skill.
3. Collects the full assistant text, extracts `datagrok-exec` fenced
   blocks, and runs each prompt's rubric:
   - `expected_path` — `datagrok-exec` block emitted?
   - `expected_helpers` — every regex matches in at least one block?
   - `expected_symbols` — every regex matches in at least one block?
   - `forbidden_symbols` — zero regexes match in any block?
4. Writes a markdown report with summary, per-prompt table, and deep
   failures.

## Test mode rationale

We call the SDK directly rather than spinning the container because:
- Same code path Claude actually executes (the container is a thin
  Hono/WebSocket wrapper around the SDK; the prompt assembly is mirrored
  in `src/prompt-assembly.ts`).
- No Docker rebuild between iterations.
- No WebSocket / MCP / user-files complexity — we want to isolate the
  *prompt* effect.
