# Healing API

> Status: full detail. This is the surface test authors see and depend on.

## Module

```ts
import * as healing from '@datagrok-libraries/self-healing-locators';
```

The library is framework-agnostic at the type level. Concrete adapters
ship for Puppeteer (used by Datagrok package tests) and Playwright
(used by the standalone E2E suite). The core API works the same way
through either adapter.

## Configuration

Done once per test process, typically in a global setup file.

```ts
healing.configure({
  registryPath: 'libraries/self-healing-locators/registry/v1',
  auditLogPath: 'test-output/self-healing/locator-events.jsonl',
  queuePath:    'test-output/self-healing/healing-queue.jsonl',
  thresholds: { high: 0.80, mid: 0.55 },
  enableScreenshots: false,
  enableTier4Visual: false,
  primarySelectorTimeoutMs: 1000,
  perTierBudgetMs: { tier1: 200, tier2: 300, tier3: 200, tier4: 500 },
  // Optional overrides; see docs/03-confidence-model.md for defaults
  weights: undefined,
});
```

`configure()` is idempotent. Calling it twice with identical options is
a no-op. Calling it with conflicting options throws — config is global
state, but it's not mutable global state.

## Core: `resolve()`

The drop-in replacement for `page.locator(selector)`.

```ts
const element = await healing.resolve(page, '[data-testid="submit"]', {
  anchorName?: string;          // explicit name; auto-derived otherwise
  action?: 'click' | 'type' | 'assert' | 'read' | 'hover';
  purpose?: 'normal' | 'destructive' | 'auth';  // policy hint
});
```

Returns the framework's element handle (Puppeteer `ElementHandle`,
Playwright `Locator`). Callers use it normally — `await element.click()`,
etc.

Behavior:

1. Tries the primary selector (with the configured timeout).
2. On miss, runs the resolver chain.
3. Returns the resolved element on HIGH or MID; throws on LOW or fail.
4. Writes one audit event per call (success or fail).

If the primary selector hits, `resolve()` is just a passive
fingerprint-refresh wrapper with negligible overhead.

## Explicit anchors: `anchor()`

Preferred for critical paths. Forces an explicit name and surfaces
metadata into the fingerprint.

```ts
const submit = await healing.anchor(page, 'login-submit', {
  selector: '[data-testid="submit"]',
  action: 'click',
  purpose: 'auth',
});
await submit.click();
```

Equivalent to `resolve()` with `anchorName: 'login-submit'`, but:

- The `anchorName` is required, not derived. Renames go through a
  deprecation path (see § "Renaming anchors").
- The `purpose` field is recorded in the fingerprint and influences
  the policy gates (e.g. `auth` requires explicit anchors; the policy
  refuses to heal auto-derived anchors on auth flows).
- The library can attach optional fingerprint hints
  (`hints.elementIsInsideViewer`) without the author having to compute
  them.

## Renaming anchors

Anchor names are part of the contract between tests and the registry.
A rename is a deliberate refactor, not a silent fix.

```ts
healing.deprecateAnchor('old-name', { until: '2026-09-01' });
const submit = await healing.anchor(page, 'new-name', { ... });
```

The library:

- Loads both fingerprints; new anchor inherits from old until the
  deadline.
- Writes a warning to the audit log on each use of the old name.
- After the deadline, the old anchor is pruned by `grok-heal registry prune`.

## Reading the audit log

For dashboards, custom reporters, or test-time introspection.

```ts
const events = await healing.audit.read({ since: '2026-04-01' });
const counters = healing.audit.snapshot();
// counters: { total_resolves, primary_hits, tier1_hits, ..., failed }
```

The audit log is jsonl, append-only. The library exposes a streaming
reader so consumers do not have to parse it themselves.

## Marking review decisions

After a healing PR is reviewed, decisions feed back to the offline tool
so it doesn't re-propose the same heal next run.

```ts
// Programmatic alternative to `grok-heal mark`
await healing.audit.markCase(caseId, {
  decision: 'accepted' | 'rejected',
  reason?: string,
});
```

Programmatic marking is mostly for the CLI; humans typically use
`grok-heal mark` from the command line.

## Error types

```ts
class HealingError extends Error {
  caseId: string;
  reason:
    | 'no_fingerprint'
    | 'all_tiers_missed'
    | 'ambiguous_candidates'
    | 'low_confidence'
    | 'platform_version_skew'
    | 'registry_corrupt'
    | 'unhealable_destructive'
    | 'unhealable_auth';
  details: { /* per-reason fields */ };
}
```

`HealingError` extends `Error` so existing test reporters render it
sensibly. The `caseId` lets a developer pull up the audit entry directly
from the failure message.

## What the API never does

- **Never returns null silently.** Either resolves to an element or
  throws. There is no in-band "didn't find it" sentinel.
- **Never retries on transient failures.** That's the test framework's
  job. Self-healing fixes locator drift, not flakiness.
- **Never logs to stdout.** All logging goes through the audit log or
  a configurable logger callback.
- **Never makes network calls.** No metrics services, no telemetry, no
  Anthropic SDK in this package.

## Adapter packages

The core library exposes interfaces; adapters implement them.

```
@datagrok-libraries/self-healing-locators                  // core
@datagrok-libraries/self-healing-locators-puppeteer        // adapter
@datagrok-libraries/self-healing-locators-playwright       // adapter
```

Adapters are thin: they translate between core types and the
framework's element/page types, and they implement the framework-specific
parts of fingerprint capture (taking a screenshot, reading the
accessibility tree, etc.).

The Datagrok Puppeteer test framework uses the Puppeteer adapter. The
Playwright E2E suite uses the Playwright adapter. The core API surface
is identical from the test author's perspective.

## API stability

- `1.0` is the first stable release; signatures above ship as is.
- Pre-1.0 versions may evolve based on early-adopter feedback from the
  volunteer package(s) in Stage 2 of the rollout.
- Breaking changes after 1.0 follow standard semver and require a
  deprecation window of at least one minor version.
