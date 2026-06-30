---
status: complete
phase: 12-spectronaut-input-coverage
source: [12-VERIFICATION.md]
started: 2026-05-15T15:17:36Z
updated: 2026-05-15T19:06:19Z
---

## Current Test

[testing complete]

## Tests

### 1. Re-import the 2.6 GB reference TSV after the malformed/filtered split (12-04)
expected: Open `~/Downloads/2026-05-13 BP DMD WT.tsv` via **Proteomics | Import | Spectronaut Report**. Import completes with no V8 string-length or OOM error; the TaskBar progress bar advances monotonically; no Chrome "Page Unresponsive" dialog appears; switching browser tabs works mid-import; the resulting DataFrame has approximately 8,328 protein rows and 24 sample columns; tags `proteomics.source=spectronaut`, `proteomics.preNormalized=true`, `proteomics.groups` with 2 conditions are set; the full pipeline Annotate → Normalize → Impute → DE → Volcano runs to completion. **New 12-04 acceptance:** NO false "skipped N malformed line(s)" message fires on this otherwise-correct import — by-design-filtered rows (CON__/REV__ decoys, numeric q > threshold) must be silently dropped, matching the text path.
result: pass — user confirmed re-import of the 2.6 GB reference file fires no false "malformed line(s)" message after the 12-04 split (approved 2026-05-15)

## Summary

total: 1
passed: 1
issues: 0
pending: 0
skipped: 0
blocked: 0

## Gaps

The OOM-prevention, monotonic-progress, tab-responsiveness, and full-pipeline aspects of this test were confirmed PASS in the pre-12-04 run (recorded 2026-05-15T16:35:12Z). The 12-04 gap closure changed the expected observable: by-design-filtered rows must no longer be reported as "malformed". That correction was made and merged AFTER the prior pass, so the no-false-malformed-signal criterion has not yet been confirmed on the real 2.6 GB file. Re-run is scoped to confirming that single corrected behavior alongside the unchanged import path.
