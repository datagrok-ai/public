# 07 — Healing Policy

> Status: full detail. This document defines the boundary between
> "drift we heal" and "real change we let fail." Reviewers, please
> scrutinize the examples — they encode our editorial judgment.

## 1. The principle

Self-healing repairs **how** a test finds an element. It does not repair
**what** the test asserts about that element. If the element is gone, or
its role has changed in a way that breaks the test's intent, the test
should fail. The system's job is to remove the noise of locator drift,
not to keep the test green at any cost.

This document is the answer to "what counts as drift, and what counts as
a real change?" We answer it concretely, with examples.

## 2. What we heal (drift)

A change qualifies as drift when **all** of the following are true:

1. The user-visible behavior the test is exercising is unchanged.
2. The element fulfilling that behavior still exists, with the same
   semantic role.
3. The change is in the locator path: id renamed, classes hashed
   differently, structural wrapper added, accessible name normalized,
   etc.
4. The fingerprint can be re-attached to a single element on the new
   page with confidence ≥ MID.

If any of those is false, we don't heal.

### Examples we heal

| Change | Why we heal |
|---|---|
| `data-testid="submit"` → `data-testid="submit-btn"` | Same role, same name, same position. Pure renaming. |
| Class `.btn-primary` replaced by `.css-1a2b3c4` | CSS-in-JS migration. Stable classes filtered out; semantic match remains. |
| Button moved into a new `<form-actions>` wrapper | Structural shift. Tier 2 finds it via the parent chain shifted by one. |
| `aria-label="Submit form"` → `aria-label="Submit"` | Fuzzy match within Levenshtein threshold; role unchanged. |
| Viewer's toolbar reordered | bounding box changed, but `dgWidgetType + dgViewerType + accessible name` still match. |
| `<button>` replaced by `<div role="button">` | Role match still works; tag mismatch is small structural penalty. |
| Two scrollable panes swap order | Tier 2 finds the matching pane by stable classes + parent chain. |

### Examples we heal **with a flag** for review

When confidence lands in the MID band, we heal but ask a human to
double-check. Examples that typically land here:

- Multiple structural changes at once (wrapper added AND class hash
  changed AND parent chain depth changed).
- The Levenshtein-fuzzy text match was the deciding factor.
- The accessible name changed from a verb to a synonym ("Submit" →
  "Send"). The role is the same; the wording shift is small but real.
- The visual hash similarity carried weight (visual layer activated
  because semantic and structural were inconclusive).

Flagging is the explicit acknowledgment that this is a judgment call
and a human should rubber-stamp it.

## 3. What we do not heal (real changes)

A change is real, and the test must fail, when **any** of the following
is true:

1. The element no longer exists on the page in any form the resolver
   can match with confidence ≥ MID.
2. The semantic role of the element changed (`button` → `link`,
   `dialog` → `popover`).
3. The Datagrok widget type changed (`Viewer` → `Dialog`,
   `scatter-plot` → `bar-chart`).
4. The element's intent changed: button label "Submit" became "Delete"
   (text similarity will not save this — semantic mismatch is too
   large, and the action mismatch gate may also fire).
5. Multiple top candidates are tied within the ambiguity threshold and
   we cannot pick safely.

### Examples we do not heal

| Change | Why we let the test fail |
|---|---|
| Submit button removed (refactor mistake) | Element absent. Fail. |
| "Submit" button became "Delete" | Same role and position, but the test's intent doesn't survive. The visible-text mismatch + the fingerprint's `purpose` (if explicitly set) drive this to LOW. |
| Login dialog became a fullscreen view | Widget type changed. Fail. |
| The "Run" button is now disabled by default | Element exists, but action mismatch (test does click; element is non-interactive) fires. Fail. |
| Two identical "Save" buttons appear (one in toolbar, one in dialog) | Ambiguity gate fires. Fail rather than guess. |

## 4. The gray zones we have to call

Not everything is clean. We make the following editorial choices and
document them so the team is aligned.

### 4.1 Synonym wording changes

"Submit" → "Send", "Apply" → "Save", "OK" → "Confirm". These pass our
fuzzy-match threshold sometimes and not others depending on the rest of
the fingerprint. Our policy: when text is the **deciding** signal, we
flag for review. We never silently heal a wording change.

### 4.2 Role normalization

Sometimes a `<div onclick>` becomes a `<button>`. The role changes from
implicit to explicit, but the user experience is unchanged or improved.
Our policy: heal at HIGH if `accessible name + position` match strongly
and the role change is a **strengthening** (less semantic to more
semantic). The reverse — `<button>` becoming `<div onclick>` — flags
for review because it suggests an accessibility regression worth
noticing.

### 4.3 Viewer rendering changes

Viewers (scatter plot, grid, bar chart) re-render frequently. A test
clicking a specific bar in a chart should not be a self-healing target
in the first place — those tests should use the viewer's API
(`viewer.dataFrame.currentRow = N`), not raw clicks. Our policy: if a
fingerprint sits **inside** a viewer at a deep level, we heal at the
viewer boundary (the toolbar, the title, the tab) but **refuse** to
heal at the rendered-data level. The fingerprint capture utility warns
when an anchor sits below the viewer rendering boundary.

### 4.4 Dynamic content

A heading "Welcome, Alice" tested from a session that runs as Alice. If
the test reruns as Bob, the text is "Welcome, Bob." Our policy: detect
dynamic patterns at capture (`hints.textIsLikelyDynamic`), halve the
text feature weights at scoring time, and rely on role + position. If
the test author wrote a fragile assertion, we don't fix that — the
assertion itself will catch the mismatch.

### 4.5 i18n

A button labeled "Submit" in en-US, "Senden" in de-DE. The test's locale
is fixed per run, so the fingerprint matches its locale. If the test
is parameterized over locales, the test author should use the
non-text features as anchors (role, testid, structural). If they
haven't, locale changes look like text drift. Our policy: detect
mismatched language between fingerprint and live page (via simple
language ID) and **fail** rather than heal — this is almost always a
test setup issue, not drift.

## 5. Hard rules (regardless of confidence)

These bypass the score entirely:

1. **No heal across widget type boundaries.** If a fingerprint says
   "button inside viewer X" and the candidate is "button inside dialog
   Y", do not heal even at HIGH score.
2. **No heal of destructive-action elements without a flag.** Anchors
   tagged `purpose: 'destructive'` (delete, drop, remove, etc.) always
   flag for review even at HIGH score. The cost of healing wrongly is
   higher.
3. **No heal of authentication-related elements without an explicit
   anchor.** Auto-derived anchors on login/auth flows are too risky.
   The library refuses to heal these unless the test author committed
   to an explicit `healing.anchor()`.
4. **No heal when running against an unfamiliar platform version.** If
   `runMetadata.platformVersion` is more than two minor versions ahead
   of the fingerprint's platform version, the resolver returns LOW
   automatically. The platform is allowed major UI overhauls between
   minor versions, and healing across them risks masking redesigns.

These rules are encoded as code, not config, because they are safety
properties — turning them off changes the system's character, not its
tuning.

## 6. What the policy implies for test authors

Authors are still responsible for:

- Writing meaningful assertions. We don't heal those.
- Choosing semantic selectors when possible. Tier 1 is most reliable.
- Marking destructive actions with `purpose: 'destructive'` so the policy
  applies.
- Using `healing.anchor()` for auth flows and other high-risk paths.

The library does not absolve authors of these responsibilities. It
absorbs the parts of test maintenance that don't add information.

## 7. What the policy implies for reviewers

When a healing PR lands in your queue:

- Each row in the table is one healed selector. The score and tier tell
  you how strong the evidence was. HIGH is "I would have done the same
  rename;" MID is "the system thought this was right but please look."
- The audit references in the PR body link to per-case detail. You can
  see the fingerprint, the proposed candidates, and which one won.
- If you disagree with a heal: revert the line, mark the case as
  `human-rejected` via `grok-heal mark`, and the case won't be
  proposed again until the fingerprint or page state changes.
- If you see patterns of bad heals: that's signal worth raising. The
  weights or thresholds may need tuning, or the platform may need
  better test IDs.

## 8. What the policy implies for the platform team

The cleanest lever the platform team has to reduce false heals is
**stable test IDs on UI primitives** (see `06-datagrok-integration.md`).
Every `data-testid` we get for free in the platform is one less case
the resolver has to reason about structurally.

The system will work without it; it will work better with it.

## 9. Open policy questions (please weigh in during review)

1. **MID-band default.** Heal-and-flag, or fail-and-flag? We propose
   heal-and-flag because the test exercising the rest of the flow is
   valuable even if the locator is uncertain. Reviewers may prefer
   fail-and-flag for higher safety.
2. **Auto-rejection memory.** When a reviewer rejects a heal, how long
   do we remember it? Until the fingerprint changes? Until the page
   snapshot changes? Forever? We propose: until fingerprint OR page
   snapshot hash changes.
3. **Bulk approvals.** If a release introduces 50 heals across 10
   packages, do we open one PR or ten? We propose ten (per package,
   per CODEOWNERS), but reviewers may prefer one mega-PR for visibility.
