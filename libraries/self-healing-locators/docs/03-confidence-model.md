# 03 — Confidence Model

> Status: full detail. The confidence formula gates every healing decision.
> Weights are calibratable; the formula's shape is not.

## 1. What the confidence score is — and isn't

It is: a number in `[0, 1]` indicating how strongly the candidate element
matches the recorded fingerprint, computed from objective feature
overlaps.

It is not: an estimate from a language model. Claude does not produce a
confidence number that we use directly. When the offline tool processes
LLM proposals, the same scorer that runs at runtime re-scores each
proposal against the live page. This is intentional — it puts the LLM
and the runtime in the same evaluative frame and prevents the LLM from
"talking itself" into a heal.

## 2. Formula

For a candidate element matched against a stored fingerprint:

```
score = sum_over_features( w_i * m_i ) / sum_over_features( w_i )

where:
  w_i = configured weight for feature i
  m_i = match value in [0, 1] for feature i (0 if feature is absent
        in either the candidate or the fingerprint)
```

This is a **weighted Jaccard-style ratio**, not an unbounded sum.
Properties we get for free:

- Score stays in `[0, 1]`. No artificial caps or normalizations later.
- Missing features on either side don't penalize disproportionately —
  they just don't contribute weight.
- Adding new features changes the denominator too, so an existing
  fingerprint with a missing new feature isn't suddenly downscored.

## 3. Default weights

Default weights, by tier. These are **starting values**. The eval suite
will calibrate them on real Datagrok cases before the system is trusted
in production.

### Tier 1 — semantic

| Feature | Weight | Match function |
|---|---|---|
| `testId` exact | 1.00 | binary |
| `domId` exact (passes stability filter) | 0.90 | binary |
| `ariaLabel` exact | 0.85 | binary |
| `ariaRole` + accessible `name` exact | 0.85 | binary |
| `dgWidgetType` + `dgViewerType` exact | 0.80 | binary |
| `ariaLabel` fuzzy (Levenshtein ≥ 0.85) | 0.55 | proportional to similarity |

### Tier 2 — structural

| Feature | Weight | Match function |
|---|---|---|
| `tagName` exact | 0.20 | binary |
| `stableClasses` Jaccard | 0.40 | size-of-intersection / size-of-union |
| `parentChain` match (any depth) | 0.35 | per-link partial credit |
| `indexInParent` match | 0.10 | binary, last resort |

### Tier 3 — text

| Feature | Weight | Match function |
|---|---|---|
| `visibleText` exact (after normalization) | 0.70 | binary |
| `visibleText` fuzzy (Levenshtein ≥ 0.85) | 0.45 | proportional |
| `placeholder` exact | 0.40 | binary |
| `title` exact | 0.30 | binary |

When `hints.textIsLikelyDynamic` is true, these weights are halved. We
do not trust dynamic text as a strong identifier.

### Tier 4 — visual / positional

| Feature | Weight | Match function |
|---|---|---|
| `boundingBox` IoU | 0.25 | Intersection-over-Union ratio |
| `visualHash` similarity | 0.20 | Hamming distance to similarity, threshold ≥ 0.9 |

When the screenshot feature flag is off, the visual layer contributes
nothing. The denominator shrinks accordingly.

## 4. Thresholds

The decision gate has three bands:

| Band | Range | Action |
|---|---|---|
| HIGH | `score ≥ 0.80` | Heal silently. Audit log, no flag. |
| MID | `0.55 ≤ score < 0.80` | Heal AND flag for human review. Queue offline. |
| LOW | `score < 0.55` | Fail the test. Queue offline. |

Additional gates that override the band (these are hard rules, not score
adjustments):

- **Ambiguity rule.** If two candidates score within 0.10 of each other
  in HIGH or MID, demote both to LOW. We never silently pick between
  near-tied candidates.
- **Semantic mismatch.** If `ariaRole` of the candidate differs from the
  fingerprint's `ariaRole` (e.g. fingerprint says `button`, candidate is
  `link`), demote the candidate by one band. This catches "found
  something that looks similar but is the wrong kind of thing."
- **Action mismatch.** If the test action is `click` and the candidate
  is not interactive (no `tabindex`, no `onclick`, not a focusable role),
  demote by one band.
- **Datagrok widget mismatch.** If `dgWidgetType` is recorded and the
  candidate's widget type differs, demote by one band. The platform's
  semantic layer is too informative to ignore.

## 5. Why these defaults

A few non-obvious choices:

**Tier 1 dominates.** If `testId` matches, the score is ≥ 1.00 / (sum of
present features), which alone is enough for HIGH on a typical
fingerprint. This is the right outcome — when a stable semantic anchor
matches, we should be highly confident regardless of structural drift.

**Structural is the smallest tier in aggregate.** Sum of structural
weights ≈ 1.05. Structure is a tiebreaker, not a primary signal. SPAs
restructure constantly.

**Text is heavier than structure but conditional.** Strong when stable;
halved when dynamic. The conditional logic prevents
`visibleText='Welcome, Alice'` from anchoring a test that runs as Bob.

**Visual is the lightest.** Visual signals are noisy (rendering,
fonts, themes, viewport sizes vary). They're a sanity check, not a
foundation.

## 6. Calibration

Defaults will be calibrated against the eval fixtures (see
[`09-eval-and-rollout.md`](./09-eval-and-rollout.md)). The calibration
target: maximize true heals at HIGH while keeping false heals at HIGH
below a configured budget (proposal: < 1%).

The eval suite reports per-tier contribution, so we can see whether,
say, `dgWidgetType` is overweighted in practice or whether
`stableClasses` is doing nothing useful and can be deprioritized.

Weights live in `self-healing.config.json`; calibration changes ship as
config commits, not code commits.

## 7. Why a formula and not a model

We considered training a small classifier or learned reranker. Rejected
for now because:

- We don't have labeled data yet.
- Explainability matters: a reviewer asking "why did this heal?" gets a
  per-feature breakdown from the formula. A learned model would need
  separate explainability tooling.
- The formula has fewer than 15 numbers to tune. A model would have many
  more, with attendant overfitting risk on a small corpus.
- If the eval suite shows the formula is fundamentally inadequate, we
  swap the scorer behind the same `Scorer` interface. The architecture
  permits this.

## 8. How the offline tool uses the same scorer

The offline CLI processes a queued case as follows (full details in
[`05-offline-flow.md`](./05-offline-flow.md)):

1. Load the fingerprint and the live page (via headless browser).
2. Ask Claude for up to 5 candidate selectors.
3. For each candidate, find the element on the page.
4. **Re-score the candidate using the runtime scorer** against the
   fingerprint.
5. Rank candidates by their runtime score, not by Claude's preference
   order.
6. Apply the same band gates as runtime.

This means a "great" Claude proposal that scores LOW under our scorer
gets rejected. A "weak-looking" Claude proposal that, on inspection,
recovers Tier 1 features the resolver missed gets accepted. The scorer
is the source of truth.

## 9. Audit log fields related to scoring

Every heal event includes:

```json
{
  "score": 0.83,
  "band": "HIGH",
  "matched_features": [
    {"name": "ariaLabel", "weight": 0.85, "match": 1.0},
    {"name": "ariaRole+name", "weight": 0.85, "match": 1.0},
    {"name": "stableClasses", "weight": 0.40, "match": 0.66}
  ],
  "absent_features": ["testId", "domId", "dgWidgetType"],
  "demotions_applied": [],
  "ambiguity_check": "passed"
}
```

Reviewers can reconstruct the decision exactly. No black box.
