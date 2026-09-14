# Codex repository guidance

The shared instructions for this repository live in [CLAUDE.md](CLAUDE.md) and
`.claude/`. Read and follow them when using Codex:

- Read each applicable nested `CLAUDE.md` along the path to the files being inspected
  or changed. More specific guidance takes precedence; explicit user instructions
  take precedence over repo guidance.
- Read `.claude/rules/code-style.md` and any other `.claude/rules/` files whose
  `paths` match the task. Rules without path restrictions apply throughout the repo.
- Check `.claude/skills/` for a matching workflow and read its `SKILL.md` and needed
  references before using it. Reuse these files rather than maintaining copies.
- Read `.claude/settings.json` and relevant referenced hooks. They do not execute
  automatically in Codex; follow their applicable constraints with available tools.
  In particular, regenerate `.g.ts` and `.api.g.ts` with `grok api` instead of editing them.
- When this checkout is the core repository's `public/` submodule, also read the
  enclosing repository's `AGENTS.md`, `CLAUDE.md`, and relevant `.claude/` guidance.
  Keep changes and Git operations scoped to the appropriate repository.

For BDD work, read [libraries/bdd/CLAUDE.md](libraries/bdd/CLAUDE.md), its README, and
the target package's `bdd/README.md`. Generated BDD specs are committed artifacts:
edit features or bindings, run `grok-bdd compile`, and check for drift before running.
