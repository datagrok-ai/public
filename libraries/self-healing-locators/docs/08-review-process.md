# 08 — Review Process

> Status: compact. The audit log and PR conventions are spelled out in
> earlier docs; this document just orients reviewers and authors.

## Two review surfaces

**1. Healing PRs** opened by `tools/self-healing-cli`. One per package,
auto-assigned to the package's CODEOWNERS. Body includes a per-case
table (test, anchor, old → new selector, score, tier), unhealed cases,
cost, and audit references.

**2. Review queue.** A separate listing of MID-band runtime heals
(healed but flagged) and offline cases the tool couldn't heal. Renders
from the merged audit log and the offline report. **Out of scope for
this PR**; design proposed here as a placeholder so we can cross-link.

## What a reviewer does on a healing PR

- Skim the table. Anything in the MID band gets a closer look.
- Click through to the audit reference for a case where the heal feels
  wrong. The audit has the fingerprint, the candidates Claude proposed,
  the runtime scorer's per-feature breakdown, and the winner.
- Approve, request changes, or reject specific lines.
- If a heal is wrong, revert the line and run
  `grok-heal mark <case-id> --rejected` to remember.

## What a reviewer should NOT do

- Do not approve the PR without skimming MID-band rows.
- Do not assume HIGH-band rows are always correct — they almost always
  are, but the audit is short and worth a glance for any visually
  surprising selector change.
- Do not merge a healing PR while a related platform release is mid-flight.
  Wait for the platform to settle, run the tool again, then merge.

## When to escalate

- Sustained increases in MID-band heals across releases → tune weights
  or push for `data-testid` adoption in the platform.
- A specific test consistently produces low-confidence heals → that
  test's anchors are likely too weak. Add explicit `healing.anchor()`
  with stronger fingerprint hints.
- Claude proposes implausible candidates repeatedly → the prompt or
  examples need work; raise an issue against the prompt files.

## Ownership

- **Library code:** owners of `libraries/self-healing-locators`.
- **CLI code:** owners of `tools/self-healing-cli`.
- **Healing PRs:** package CODEOWNERS by default, configurable.
- **Prompt files:** the same owners as the library, with QA team as
  required reviewer (so prompt changes get more eyes than code changes
  on the library internals).
- **Confidence weights and thresholds:** library owners + QA team
  required reviewer.
