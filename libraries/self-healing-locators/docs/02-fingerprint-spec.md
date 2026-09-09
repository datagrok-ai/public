# 02 — Fingerprint Specification

> Status: full detail. The fingerprint schema is the contract between
> capture, runtime resolution, and the offline LLM. Lock this carefully.

## 1. What a fingerprint is

A fingerprint is a structured snapshot of an element captured when a test
last passed. It contains enough orthogonal signals that, when the primary
selector breaks, we can re-identify the same element through a different
combination of signals — and quantify our confidence in the match.

A fingerprint is not a CSS selector. It is the **evidence** we use to find
or generate a working selector.

## 2. Schema (informal)

The formal JSON Schema lives in
[`schemas/fingerprint.schema.json`](../schemas/fingerprint.schema.json).
Here is the human-readable form, with rationale per field.

```ts
interface ElementFingerprint {
  // Identity
  schemaVersion: '1';
  anchorName: string;            // logical name from healing.anchor()
                                 // or auto-derived from selector + test
  capturedAt: string;            // ISO 8601
  capturedBy: 'auto' | 'explicit';
  testRef: {
    file: string;                // path relative to repo root
    suite?: string;
    test: string;
    selectorAtCapture: string;   // the selector that worked
  };

  // Tier 1 — semantic
  semantic: {
    testId?: string;             // data-testid value, if any
    domId?: string;              // id, only if it passes stability filter
    ariaLabel?: string;
    ariaRole?: string;           // button, dialog, tab, ...
    name?: string;               // accessible name (from a11y tree)
    dgWidgetType?: string;       // from getWidgetStatus(), e.g. 'Viewer'
    dgViewerType?: string;       // e.g. 'scatter-plot', 'grid'
  };

  // Tier 2 — stable structural
  structural: {
    tagName: string;
    stableClasses: string[];     // post-filter (see § 4)
    parentChain: ParentLink[];   // up to 3 levels, see below
    indexInParent?: number;      // last-resort sibling index
  };

  // Tier 3 — text content
  text: {
    visibleText?: string;        // normalized, trimmed, length-capped
    placeholder?: string;
    title?: string;
    valueAtCapture?: string;     // for inputs, only if not sensitive
  };

  // Tier 4 — visual / positional
  visual: {
    boundingBox?: {              // viewport-relative at capture
      x: number; y: number;
      width: number; height: number;
    };
    viewportSize?: { w: number; h: number };
    visualHash?: string;         // pHash; only when feature flag enabled
    screenshotRef?: string;      // file path under registry, optional
  };

  // Stability hints discovered at capture time
  hints: {
    domIdLooksStable: boolean;   // see filter logic § 4
    classesAreCssInJs: boolean;
    textIsLikelyDynamic: boolean; // detected via heuristics, e.g. "Hello, $user"
    elementIsInsideViewer?: string; // viewer type if applicable
  };

  // For the LLM and for debugging
  contextSnippet: {
    accessibilityNode?: object;  // a11y subtree at this element, depth 2
    htmlSnippet?: string;        // sanitized outerHTML, truncated
  };
}

interface ParentLink {
  tagName: string;
  testId?: string;
  ariaRole?: string;
  dgWidgetType?: string;
  stableClasses?: string[];      // already filtered
}
```

## 3. Why these fields and not others

Every field justifies its weight in the scoring formula. We omit:

- **XPath of the element.** Too brittle, encodes everything that drifts.
  We reconstruct an XPath at runtime if needed; we never store one as
  ground truth.
- **Full outer HTML.** Heavy, leaky, encourages bad matching. We keep a
  truncated snippet only as context for Claude in the offline phase.
- **Computed styles.** Volatile across themes and platform versions.
- **Event listeners.** Not observable reliably in headless contexts.
- **Embeddings.** Not yet. Adds infrastructure and explainability cost
  without a clear win at our current failure rates. Revisit after the
  eval suite has data.

We include:

- **`name` from the accessibility tree** even when it duplicates `ariaLabel`
  or `visibleText`. The a11y tree resolves these per the platform's
  algorithm, which is what assistive tech and Playwright `getByRole` use.
  Storing it makes Tier 1 matches more robust.
- **`dgWidgetType` and `dgViewerType`** because Datagrok has a layer of
  semantics above the DOM that most SPAs lack. Using it gives us a stable
  anchor that survives many DOM-level changes.
- **`hints`** because they are computed once at capture and used by both
  the runtime resolver and the offline prompt. Storing them avoids
  recomputation and ensures the LLM and the resolver see the same view.

## 4. Normalization rules

At capture, raw values pass through normalizers. Both capture and
resolution use the same normalizers — the registry stores the **post-
normalization** form. Normalizers are part of the library's public API
because tests may need to reproduce them.

### 4.1 `domId` stability filter

A `domId` is recorded only if it passes all of:

- length between 2 and 64
- contains at least one lowercase letter
- does not match a UUID-like pattern (`/^[0-9a-f-]{12,}$/i`)
- does not match a hash-like pattern (`/^[a-zA-Z0-9_-]{6}$/` with high
  entropy by Shannon estimate ≥ 4.5 bits/char)
- is unique on the page at capture time

If filtered, `domId` is `undefined` and `hints.domIdLooksStable = false`.

### 4.2 Class filter

A class survives into `stableClasses` only if:

- it is not in the configured CSS-in-JS pattern list (default:
  `/^css-[a-z0-9]{4,}$/`, `/^_[a-zA-Z0-9_]{6,}$/`,
  Material/Emotion-style hash patterns)
- it is not transient: not in `is-active`, `is-hover`, `is-focused`,
  `is-selected`, etc. (configurable list)
- it appears on at most 100 elements at capture time (uniqueness signal;
  configurable threshold)

The full filtered set is preserved in `hints.classesAreCssInJs` so we
know whether class-based matching is even worth attempting.

### 4.3 Text normalization

- Trim and collapse whitespace.
- Truncate to 200 characters; original length stored separately.
- Detect dynamic patterns (interpolation markers, dates, counts,
  user-specific strings) and set `hints.textIsLikelyDynamic`. The
  detector is heuristic — see `prompts/examples/` for cases.

### 4.4 Parent chain

We store **at most 3 levels up**, stopping at:

- the test view root (Datagrok view boundary)
- a `<body>` or document root
- a parent we can already identify uniquely by `testId` or `dgWidgetType`

Storing more is wasteful; storing fewer makes structural matching too
weak in deeply nested viewers.

## 5. Visual layer (feature-flagged)

By default the visual layer is **off**. Reasons:

- Doubles capture time (screenshot + hash) on every green run.
- Increases registry size by 1–10 KB per element when screenshots are
  retained, plus the visual hash itself.
- Per-deployment privacy considerations may apply to screenshot content.

When the `screenshots` feature flag is on:

- `visual.boundingBox` and `visual.viewportSize` are always captured
  (cheap, no image content).
- `visual.visualHash` is a perceptual hash (default: pHash) of the
  element-only screenshot, computed locally. No image data leaves the
  test machine if screenshots themselves are not retained.
- `visual.screenshotRef` is set only when image retention is enabled,
  pointing to a file under `registry/screenshots/<hash>.webp`.

The offline prompt **may** include a screenshot as an `image` content
block when the flag is on; the prompt template handles the absence
gracefully.

## 6. Storage

The registry is a **versioned directory of JSON files**. Recommended layout:

```
libraries/self-healing-locators/registry/
  v1/
    packages/
      <PackageName>/
        <test-file-id>.fingerprints.json   # all anchors for one test file
    e2e/
      <suite-name>/
        <test-id>.fingerprints.json
    screenshots/                            # only if flag on
      <hash>.webp
  config/
    self-healing.config.json
    css-in-js-patterns.json
```

### 6.1 Why a central registry, not co-located files

We considered three options:

- **Option A — co-located** (`my.spec.ts` + `my.fingerprints.json` next to it).
  Pros: discoverability, atomic moves with tests.
  Cons: every test directory grows; harder to scan all fingerprints across
  the repo for migrations or audits; spreads test concerns into product
  package directories.

- **Option B — central registry** (this proposal).
  Pros: clear ownership boundary; bulk operations possible; schema
  migrations live in one place; offline tools have a single source.
  Cons: rename-a-test requires updating two places; PRs that touch many
  tests touch the registry.

- **Option C — Datagrok entity** (store in the platform itself).
  Pros: collaboration features, sharing, history.
  Cons: tests need to talk to the platform to read fingerprints, breaking
  P5 (offline-capable). Also creates a circular dependency: platform tests
  depend on the platform being up.

We choose **B**. The registry is a sibling under
`libraries/self-healing-locators/registry/` so the library and its data
ship together. This is one of the open questions in the PR description —
reviewers may push for co-location, in which case we adapt.

### 6.2 File granularity

One JSON file per **test file**, not per test or per anchor. This balances:

- not too coarse (one giant file per package → merge conflicts)
- not too fine (one file per anchor → thousands of tiny files)

Within the file, anchors are keyed by `anchorName`. The capture API
guarantees stable names; see `api/healing-api.md`.

### 6.3 Versioning

Schema version is in every file. The library's loader checks the version
on read. Mismatch behavior:

- Older version, known migration → migrate in memory, log a warning, write
  back on next capture.
- Older version, no migration → fail closed. Manual intervention required.
- Newer version than library knows → fail closed. Library is too old.

## 7. Capture lifecycle

Two capture modes coexist; both produce identically-shaped fingerprints.

### 7.1 Passive (auto-capture on green runs)

When a test calls `healing.resolve(selector)` and the primary selector
**works** (no fallback needed), the resolver records or refreshes the
fingerprint for that anchor. Cost: a few milliseconds per call.

Auto-derived `anchorName`: `<testFile>::<selector>`. Stable across runs
of the same test, brittle if the test author changes the selector string
(this is fine — that's a deliberate change).

### 7.2 Explicit (`healing.anchor()`)

```ts
const submit = await healing.anchor('login-submit', page.locator('#submit'));
await submit.click();
```

The author commits to a name. Renaming is a deliberate refactor with a
deprecation path (see `api/healing-api.md`).

Use explicit anchors for:

- Critical-path elements where heal decisions need extra scrutiny.
- Elements whose default auto-derived name would be ambiguous (selectors
  built dynamically).
- Elements where the author wants to attach metadata (`role`, `purpose`)
  that informs the offline prompt.

### 7.3 What we never do at capture

- Capture during a flaky run. The library refuses to update a fingerprint
  if the runtime resolver had to fall back to any tier beyond Tier 0
  (primary selector worked exactly).
- Capture sensitive content. The text normalizer has a configurable
  redaction list; values matching it are stored as `<redacted>` and the
  fingerprint records `text.valueAtCapture` only if the field is
  whitelisted (e.g. button labels, never password fields).

## 8. Lifecycle: when fingerprints are deleted

- When a test is deleted (CI hook prunes orphaned anchors, opens
  cleanup PR).
- When `anchorName` changes (old entry kept for one minor version with
  `deprecated: true` flag, then pruned).
- Manually, via `grok-heal registry prune` for stale or known-broken
  entries.

The registry is not write-only. It must be maintainable.

## 9. Open questions for review

1. **Storage location** — central under the library, or co-located with
   tests? See § 6.1.
2. **Per-anchor TTL** — should fingerprints captured before a known
   platform breaking-change release auto-expire?
3. **Sensitive value handling** — is the redaction list approach enough,
   or do we need a per-test opt-in for `valueAtCapture`?
4. **Screenshot retention** — if the flag is on, do we store images in
   the registry (large, slow) or only ephemerally during one CI run?
