# grok-heal — Self-Healing Maintenance CLI

Offline companion to `@datagrok-libraries/self-healing-locators`. Reads
the healing queue produced by failed runtime resolutions, asks Claude
to propose new selectors for cases the deterministic resolver couldn't
handle, validates each proposal against a live page, and opens a PR
with codemod'd test source updates.

## Status

**Design phase.** This directory currently contains design docs only.
See [`docs/cli-overview.md`](./docs/cli-overview.md) and the library's
[offline flow doc](../../libraries/self-healing-locators/docs/05-offline-flow.md)
for the proposed behavior.

## Commands (planned)

```
grok-heal queue              List cases waiting for offline processing.
grok-heal run [--model auto] Process the queue end to end.
grok-heal validate <case-id> Dry-run a single case (no PR).
grok-heal pr                 Open or update healing PR(s).
grok-heal mark <case-id>     Record a review decision on a case.
grok-heal eval               Run the regression eval suite.
grok-heal registry prune     Clean up stale fingerprints.
```

## Configuration

Per-repository config lives at
`tools/self-healing-cli/config.json` and may be overridden by env vars:

```json
{
  "model_selection": {
    "default": "claude-haiku-4-5",
    "structural": "claude-sonnet-4-6",
    "escalation": "claude-opus-4-7"
  },
  "cost_cap_usd_per_run": 10.0,
  "case_retry_limit": 2,
  "cache_ttl_days": 7,
  "pr_grouping": "per-package",
  "pr_default_reviewer": "CODEOWNERS",
  "browser_engine": "playwright",
  "headless": true,
  "registry_path": "../../libraries/self-healing-locators/registry/v1",
  "queue_path": "test-output/self-healing/healing-queue.jsonl"
}
```

The Anthropic API key is read from `ANTHROPIC_API_KEY`. The CLI fails
loudly if it isn't set; it never falls back to a different provider.

## Out-of-scope (and intentionally so)

- Mutating the fingerprint registry. Fingerprint updates are a
  capture-time concern handled by the runtime library.
- Modifying tests outside of the selector argument at the recorded
  call site.
- Auto-merging PRs.
- Calling any LLM other than Claude.

## License

Same as the parent repository.
