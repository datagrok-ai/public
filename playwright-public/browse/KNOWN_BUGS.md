# Browse autotests — known platform bugs

This document tracks platform bugs that are **knowingly reproduced** by Browse autotests via
the `test.fail()` Playwright annotation. The tests run end-to-end, the platform throws the
documented error, and Playwright records the failure as **expected** — counted as ✓ passed.

When a bug is fixed on the server, the corresponding test starts producing an **unexpected
pass**. Playwright then reports `test failed (unexpected status)`, which signals to the
maintainer that the `test.fail()` annotation should be removed (and the manual case marked
as resolved).

## File

`playwright-tests/e2e/browse/modelhub.test.ts`

## Honest (non-suppressed) reproductions

These tests assert clean behavior and go **red** whenever the bug fires — the error is *not*
suppressed. They are listed here so reviewers can attribute the failure to a known ticket.

### Browse-DemoApps — Bioinformatics / Similarity, Diversity

- **File:** `playwright-tests/e2e/browse/demo_apps.test.ts`
- **Jira:** [GROK-18050](https://reddata.atlassian.net/browse/GROK-18050) — *"Demo: Bio:
  Similarity, Diversity: random clicking caused 'NullError: method not found: 'gdi' on null'"*.
  The ticket is marked **Done**, but the bug **reproduces again on dev** — a regression.
- **How the test catches it:** after opening the demo, `pokeView()` performs the "random
  clicking" (click view body → click a grid cell → hover the viewer). One of those clicks lands
  on the broken state and Dart throws `NullError: method not found: 'gdV' on null` (the minified
  accessor name varies per build; `gdi`/`gdV` are the same method), captured by the `console` /
  `pageerror` listener in `watchErrors` → `expectNoErrors` fails.
- **Intermittent:** the click uses fixed coordinates, so it only sometimes hits the broken
  element — the test is **flaky by nature** (it mirrors the bug's own intermittence). We keep it
  as a normal assertion on purpose: a real, truthful failure. The test is tagged with the issue
  via `test.info().annotations` so the ticket shows up in the report.
- **When to remove the annotation:** once GROK-18050 is genuinely fixed on dev and the test
  stops failing across repeated runs, drop the `ref` from the `Similarity, Diversity` entry.

## Tests using `test.fail` (suppressed — counted as expected failures)

### 1. Browse-ModelHub-02 — single click on a model in the tree

- **Jira:** [GROK-19740](https://reddata.atlassian.net/browse/GROK-19740)
- **Expectation per manual case:** clicking a model in `Apps > Compute > Model Hub > Uncategorized`
  selects it; Context Panel updates; no errors.
- **Current behavior on dev:** the first click triggers an exception in Dart-side code:
  ```
  TypeError: p.append is not a function
      packages/$sdk/lib/async/future_impl.dart 142:14    _FutureListener.handleError
      packages/$sdk/lib/async/future_impl.dart 648:38    _Future._propagateToListeners.handleError
      packages/$sdk/lib/async/future_impl.dart 486:5     _Future._completeError
      packages/$sdk/lib/async/future_impl.dart 55:5      _SyncCompleter._completeError
      ...
  ```
- **Annotation in code:**
  ```ts
  test.fail(true, 'Platform regression GROK-19740 — model click throws inside Dart code.');
  ```

### 2. Browse-ModelHub-03 — double click on a model

- **Jira:** [GROK-19965](https://reddata.atlassian.net/browse/GROK-19965)
- **Expectation per manual case:** double-click on a model opens it as a persistent view;
  right-click → Run remains a working fallback.
- **Current behavior on dev:** identical stack to ModelHub-02 — same `TypeError: p.append is
  not a function` from `_FutureListener.handleError`. The double-click code path lands on the
  same broken click handler.
- **Annotation in code:**
  ```ts
  test.fail(true, 'Platform regression GROK-19965 — double-click on a model throws in Dart code.');
  ```

## How `test.fail` works in Playwright

| Test outcome              | What Playwright reports |
| ------------------------- | ----------------------- |
| Test throws (as expected) | ✓ counted as **passed** (expected failure)             |
| Test passes (bug fixed!)  | ✗ counted as **failed** (unexpected pass) — fix needed |

## When to remove the annotation

1. The team merges the fix for GROK-19740 / GROK-19965 on dev.
2. The next CI run of the Browse suite shows `unexpected pass` for one of the ModelHub tests.
3. Remove the `test.fail(...)` line from that test in `modelhub.test.ts`. The test then runs
   as a normal regression assertion.
4. Update the corresponding manual case in `browse_manual_tests2.md` (mark ref as resolved).

## Related context

- These tests intentionally stay in the suite (rather than being skipped) so that:
  - The regression is automatically picked up the moment the platform fix lands.
  - Reviewers / CI logs continuously surface the known-bug attribution.
- For the same reason, the manual case `Browse-ModelHub-02` / `Browse-ModelHub-03`
  expectations are *not* relaxed; we want them to assert clean behavior once the platform
  catches up.
