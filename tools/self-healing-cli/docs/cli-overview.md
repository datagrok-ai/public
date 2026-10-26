# CLI Overview

> Status: compact. Detailed pipeline lives in
> `libraries/self-healing-locators/docs/05-offline-flow.md`.

## What `grok-heal` is

The offline maintenance binary for the self-healing system. It consumes
the healing queue, calls Claude on cases the deterministic runtime
resolver couldn't heal, validates proposals against a live page, and
opens a PR with the resulting selector changes.

## Where it runs

- **CI (scheduled).** Nightly job after the day's test runs settle.
  Best for steady-state maintenance.
- **CI (post-release).** Triggered after a platform release that
  produces queued cases.
- **Locally.** Engineers can run a single case through the pipeline
  during investigation: `grok-heal validate <case-id>`.

## Commands

### `grok-heal queue [--filter ...]`

Prints pending cases. Filters: `--package`, `--anchor`, `--since`,
`--band` (mid|low|fail). No side effects.

### `grok-heal run`

End-to-end processing of the queue. Honors the cost cap. Produces a
report and (unless `--dry-run`) opens PR(s).

Flags:
- `--dry-run` — skip codemod and PR opening.
- `--cost-cap-usd <N>` — override config.
- `--model auto|haiku|sonnet|opus` — force a specific model. Default `auto`
  uses the rules in the offline-flow doc § 5.
- `--package <name>` — restrict to one package's cases.

### `grok-heal validate <case-id>`

Reproduces one case end-to-end. Prints the Claude proposals, the
validator's per-candidate scores, the chosen winner if any, and the
diff that would be applied. Never mutates anything. Used during
investigation and during eval calibration.

### `grok-heal pr [--package <name>]`

Opens or updates the healing PR(s). Idempotent: if a PR already exists
for the package and the new content is unchanged, no action.

### `grok-heal mark <case-id> --accepted|--rejected [--reason <text>]`

Records a review decision. `--rejected` cases are remembered until
either the fingerprint or the page snapshot hash changes; the same
heal will not be re-proposed in the meantime.

### `grok-heal eval [--fixtures <path>]`

Runs the regression eval. Used in CI to gate changes to the library,
the prompt, the weights, or the model selection rules. Reports per-tier
contribution and per-band metrics; see
`libraries/self-healing-locators/docs/09-eval-and-rollout.md`.

### `grok-heal registry prune`

Removes orphaned fingerprints (anchors no longer used by any test) and
deprecated anchors past their deadline. Opens a cleanup PR.

## Outputs

Per run:

- `tools/self-healing-cli/output/healing-report-<run-id>.json` —
  summary, costs, per-case status.
- `tools/self-healing-cli/output/healing-results/<case-id>.json` —
  Claude proposals + validator scores + decision.
- One PR per package (or one combined PR per `pr_grouping` config),
  with the table of heals in the body.

## Cost control

- Hard per-run cap from config; checked before each Claude call.
- Per-case retry cap of 2 (Haiku → Sonnet → stop).
- Response cache keyed on `(fingerprint hash, page snapshot hash, model)`,
  TTL 7 days by default.
- Cases over the cap are deferred to the next run.

## Failure modes

The CLI exits non-zero only on infrastructure failures (cannot read
queue, cannot write report, cannot reach git remote). Per-case failures
land in the report with structured reasons; a queue with all cases
failing produces a `0` exit and a report describing why. CI dashboards
read the report, not the exit code.

## Dependencies

Runtime library (`@datagrok-libraries/self-healing-locators`) — for
schemas, scorer, registry I/O.

Anthropic SDK (`@anthropic-ai/sdk`) — Claude client. The only LLM
provider the CLI supports.

`playwright` (or `puppeteer`, depending on which adapter the package
uses for tests) — for candidate validation against a live page.

`jscodeshift` or equivalent — for the codemod step.

`@octokit/*` (or platform-equivalent) — for opening PRs.
