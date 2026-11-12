# 04 — Runtime Flow

> Status: full detail. This describes what happens during a single test
> run, step by step, with timing and failure modes.

## 1. Entry point

A test calls `healing.resolve(selector, options?)` instead of, or as a
wrapper around, `page.locator(selector)`. The resolve function is the
only public entry point that triggers the resolver chain.

```ts
// Drop-in for Playwright
const submit = await healing.resolve(page, '[data-testid="submit"]', {
  anchorName: 'login-submit',  // optional; auto-derived if absent
  action: 'click',             // hint for the action-mismatch gate
});
await submit.click();
```

The `healing.anchor()` helper (see `api/healing-api.md`) is a thinner
wrapper that always enforces an explicit name and is preferred for
critical-path elements.

## 2. Step-by-step flow

```
┌───────────────────────────────────────────────────────────────────┐
│ 0. Try the primary selector with a short timeout (default 1s).    │
│                                                                   │
│    Hit  → record a passive fingerprint refresh, return element.   │
│    Miss → enter the resolver.                                     │
└───────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌───────────────────────────────────────────────────────────────────┐
│ 1. Load fingerprint for `anchorName` from the registry.           │
│                                                                   │
│    Found     → continue.                                          │
│    Not found → fail loudly. Self-healing requires a fingerprint;  │
│                  we never invent one from a missing selector.     │
└───────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌───────────────────────────────────────────────────────────────────┐
│ 2. Run Tier 1 (semantic).                                         │
│                                                                   │
│    Build candidate selectors from the fingerprint:                │
│      - [data-testid="<testId>"]  (if recorded)                    │
│      - #<domId>                  (if hint says stable)            │
│      - [role="<role>"][aria-label="<label>"]                      │
│      - getByRole(role, { name })  via a11y tree                   │
│      - dg-widget-type lookup via getWidgetStatus()                │
│                                                                   │
│    For each that resolves to ≥1 element, evaluate against the     │
│    full fingerprint with the scorer. Keep all results.            │
└───────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌───────────────────────────────────────────────────────────────────┐
│ 3. Run Tier 2 (structural) only if Tier 1 best score < HIGH.      │
│                                                                   │
│    Build candidates from:                                         │
│      - tag + stableClasses intersection on the page               │
│      - parent-chain anchors found in Tier 1 (their descendants)   │
│      - sibling index within identified parent                     │
└───────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌───────────────────────────────────────────────────────────────────┐
│ 4. Run Tier 3 (text) only if best score still < HIGH.             │
│                                                                   │
│    Build candidates from:                                         │
│      - getByText(visibleText)                                     │
│      - getByPlaceholder(placeholder)                              │
│      - getByTitle(title)                                          │
│                                                                   │
│    Skip text matching when hints.textIsLikelyDynamic.             │
└───────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌───────────────────────────────────────────────────────────────────┐
│ 5. Run Tier 4 (visual) only if visual flag enabled and best score │
│    still < HIGH.                                                  │
│                                                                   │
│    Build candidates from:                                         │
│      - elements whose bounding box has IoU ≥ 0.5 with stored box  │
│      - elements whose visualHash is within Hamming threshold      │
└───────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌───────────────────────────────────────────────────────────────────┐
│ 6. Apply hard gates:                                              │
│      - ambiguity (top two within 0.10) → demote                   │
│      - semantic mismatch (role differs) → demote                  │
│      - action mismatch (non-interactive for click) → demote       │
│      - dgWidgetType mismatch → demote                             │
└───────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌───────────────────────────────────────────────────────────────────┐
│ 7. Decide:                                                        │
│      score ≥ HIGH → heal silently, return element                 │
│      MID ≤ score  → heal, return element, write 'review' flag    │
│                       to audit log AND queue for offline          │
│      score < MID  → throw, queue for offline                      │
└───────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌───────────────────────────────────────────────────────────────────┐
│ 8. Write to audit log (jsonl). Always.                            │
└───────────────────────────────────────────────────────────────────┘
```

## 3. Why short-circuit between tiers

We could compute every tier's contribution always, then sum. Instead we
short-circuit on HIGH. Reasons:

- Performance. Tier 4 (visual hashing) is by far the slowest. Skipping
  when not needed keeps the median runtime cost in the sub-100ms range
  for healed cases.
- Tier 1 hits should not need cross-tier corroboration. If `testId`
  matches and `ariaRole` matches, we don't need to also ask about the
  bounding box.

For LOW-confidence runs (score < HIGH after Tier 1), we proceed through
all enabled tiers because we genuinely need more evidence.

## 4. Performance budget

Default budget per `resolve()` call:

| Phase | Budget |
|---|---|
| Primary selector attempt | 1000 ms (configurable) |
| Tier 1 (semantic) | 200 ms |
| Tier 2 (structural) | 300 ms |
| Tier 3 (text) | 200 ms |
| Tier 4 (visual, if on) | 500 ms |
| Scorer + decision | 50 ms |
| **Total ceiling** | **~2.25 s** when all tiers run |

If a resolver phase exceeds its budget, it returns its best partial
result and continues. We never block a test indefinitely; tests time out
faster than humans do.

For the green-path (primary selector works), overhead is one fingerprint
write to memory + a deferred flush — under 5 ms typical.

## 5. The audit log

Path is configurable; default
`<test-output-dir>/self-healing/locator-events.jsonl`. One JSON line per
event. Schema in `schemas/healing-event.schema.json`. Fields:

- `eventId` (UUID)
- `timestamp`
- `testRef` (file, suite, test name)
- `anchorName`
- `originalSelector`
- `resolvedSelector` (the selector that actually worked, if any)
- `tierUsed` (1–4 or `primary` or `failed`)
- `score`, `band`
- `matched_features`, `absent_features`, `demotions_applied`
- `ambiguity_check`
- `flagForReview` (boolean)
- `queuedOffline` (boolean)
- `error` (only on failure)

The audit log is **append-only** during the run. CI uploads it as a
build artifact and merges it into the central healing queue if any
events have `queuedOffline: true`.

## 6. The healing queue

Path: `<central-location>/self-healing/healing-queue.jsonl`. Same fields
as the audit log, plus:

- `domSnapshotPath` — path to a sanitized DOM snapshot taken at the moment
  of the failure. Used by the offline tool to reproduce the page state.
- `accessibilitySnapshotPath` — path to the accessibility tree snapshot.
- `screenshotPath` — only if the visual flag is on.
- `runMetadata` — branch, commit, CI run ID, browser, viewport.

Snapshots are written next to the queue entry. They're necessary because
the offline tool runs **after** the test, against a fresh page that may
already have moved on. We need the page state at the moment of failure.

Snapshot retention: configurable, default 30 days, then pruned. Pruning
also closes any queue entries that have not been resolved offline within
the retention window — they are escalated to a human dashboard with
status `expired`.

## 7. Concurrency and parallel runs

Datagrok packages are tested in parallel (CI runs many packages
simultaneously). Several concerns:

- **Audit log writes.** Each test process writes its own `.jsonl` file
  in its own output directory. CI aggregates after the run. No locking.
- **Fingerprint registry reads.** Read-only at runtime. Safe.
- **Fingerprint registry writes.** Passive auto-capture in green runs
  could race. We use atomic file replacement (write to `.tmp` then
  rename) and last-writer-wins. Fingerprint capture is idempotent for
  unchanged elements, so collisions are usually no-ops.
- **Healing queue writes.** Each test process appends to its own file;
  CI merges. Even if an unrelated test fails the same anchor at the
  same time, both events survive — the offline tool dedupes.

## 8. Failure modes and behaviors

| Situation | Behavior |
|---|---|
| Primary selector works, fingerprint missing | Capture passively. Test passes. |
| Primary selector works, fingerprint stale | Refresh. Test passes. |
| Primary selector fails, fingerprint missing | **Fail the test loudly.** Surface a clear error: "no fingerprint to heal from." This is a setup error, not a drift case. |
| All tiers run, all miss | Fail. Queue with `tierUsed: failed`. |
| All tiers run, ambiguity gate fires | Fail. Queue with reason `ambiguous`. |
| Resolver itself throws | Fail with the original selector error wrapped. Never silently swallow. |
| Registry file corrupt | Fail with a setup error. Do not attempt heuristic recovery. |
| Audit log path not writable | Log a warning; do not fail the test. The audit log is best-effort from the test's perspective; missing audit data is a CI infrastructure problem. |

## 9. What the resolver does **not** do at runtime

- Does not call any LLM.
- Does not write to the registry on a non-green run.
- Does not retry the test action. If the heal succeeds and the click
  later fails for another reason, that's a real test failure.
- Does not modify test source.
- Does not skip assertions. If a heal returns the wrong element and an
  assertion catches it, the assertion failure stands; the audit log
  shows the heal that led to it. This is exactly what we want — a way
  to investigate when self-healing made a mistake.

## 10. Observability

In addition to the per-event audit log, the runtime exposes counters
through `healing.audit.snapshot()`:

- `total_resolves`
- `primary_hits` / `tier1_hits` / `tier2_hits` / `tier3_hits` / `tier4_hits`
- `flagged_for_review`
- `queued_offline`
- `failed`

CI publishes these to the standard test reporter. Sustained increases in
`tier3_hits` or `flagged_for_review` are signals worth investigating —
they mean the platform is drifting in ways our deterministic tiers
struggle with, which may justify a `data-testid` campaign or a config
tune.
