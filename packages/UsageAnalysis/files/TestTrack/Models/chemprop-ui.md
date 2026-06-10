---
feature: models
sub_features_covered: []
target_layer: manual-only
coverage_type: regression
pyramid_layer: integration
produced_from: extracted
original_path: public/packages/UsageAnalysis/files/TestTrack/Models/chemprop.md
extraction_date: 2026-06-10
unresolved_ambiguities:
  - block-3-container-stop-run-heavy-side-effect-on-concurrent-tests
---

# Chemprop Docker container lifecycle

Manages the **chem-chemprop** Docker container via
**Browse > Platform > Dockers**. Extracted from `chemprop.md` Block 3
because the Stop → Run cycle is a heavy-side-effect operation that
disrupts concurrent Chemprop training/apply for the duration of the
restart (~5–10 min) and is unsuitable for automated parallel CI runs.

## Block 3: Manage the chem-chemprop Docker container

1. Go to **Browse > Platform > Dockers** and locate the container
   named **chem-chemprop**.
   **Verify:** the container card surfaces with its current
   status (Running or Stopped).
2. Right-click the **chem-chemprop** container and select
   **Stop**.
   **Verify:** the container transitions to Stopped (the card
   status updates; no error surface).
3. Once the container has stopped, right-click it again and select
   **Run** to restart it.
   **Verify:** the container transitions back to Running; the
   restart completes within a bounded wait; downstream Chemprop
   training / apply continue to function (sanity-check by
   navigating back to **Browse > Platform > Predictive models**
   and observing that the saved `test_chemprop` entity is still
   present).

## Notes

- **Docker container side effects.** The `chem-chemprop` Stop+Run
  cycle affects any concurrent test against the same Datagrok
  instance that exercises the Chemprop engine (training or apply
  loses inference for the duration of the restart, ~5–10 min per
  the decision-log entry for the Chem chemprop conversion
  2026-05-12). Gate this scenario behind a single-test-instance
  precondition (skip on shared CI) or route the Stop/Run through
  the `dapi.docker.dockerContainers` programmatic API so the side
  effect is bounded and detectable.
- **Container name.** On dev the container card may appear as
  `chem-chemprop` or `chem-chemprop_1` (numeric suffix). Use a
  prefix/substring match.
