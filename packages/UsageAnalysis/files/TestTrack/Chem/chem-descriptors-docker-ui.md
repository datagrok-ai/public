---
feature: chem
realizes_atlas:
  - chem.int.descriptors-docker-bounded
realizes:
  - chem.calculate.descriptors
priority: p1
target_layer: manual-only
coverage_type: regression
companion_to: chem-descriptors-docker.md
manual_only_reason: |
  Split from chem-descriptors-docker.md (Scenario 2). The container-unavailable
  bounded-error path (GROK-17621 regression guard) CANNOT be automated on the dev
  validation environment: the chem-chem Docker container is provisioned and
  reachable there, and the automated Playwright harness cannot stop a
  platform-managed container to establish the "unavailable" precondition.

  Live MCP recon on dev.datagrok.ai 2026-08-19 confirms it: opening
  System:DemoFiles/chem/smiles.csv (1000 rows, molecule column canonical_smiles),
  Chem > Calculate > Descriptors..., selecting MolWt, and clicking OK SUCCEEDED in
  ~4.5s — MolWt column appended, column count 20 -> 25, NO error balloon. So on dev
  the error-notification and "column count unchanged" observables of Step 4 are not
  producible in a spec (the call succeeds instead of failing).

  Manual execution: run in an environment where the chem-chem container is stopped
  or absent (stop it via the Datagrok Docker management panel, or use a stand where
  it was never provisioned). Then open System:DemoFiles/chem/smiles.csv,
  Chem > Calculate > Descriptors..., select any descriptor, click OK, and confirm
  within 60 seconds a visible error balloon/dialog appears (e.g. "Timeout of
  waiting for container status change"), the UI does not hang on an indefinite
  spinner, and no descriptor columns were appended (column count unchanged).
related_bugs:
  - id: GROK-17621
    status: fixed
---

# Chem — Descriptors via Docker: container-unavailable bounded error (manual)

Manual companion to `chem-descriptors-docker.md`. Covers the container-unavailable
bounded-error contract for `Chem > Calculate > Descriptors` — the GROK-17621
regression guard — which requires a stopped/absent chem-chem Docker container that
the automated dev harness cannot establish.

## Setup

1. Use an environment where the `chem-chem` Docker container is stopped or was never
   provisioned. Stop it via the Datagrok Docker management panel if it is running.
2. Confirm the Chem package is loaded and the **Chem > Calculate** top-menu is
   registered.

## Scenario: Container-unavailable bounded error (GROK-17621 regression guard)

Steps:

1. Confirm the `chem-chem` container is not running.
2. With `System:DemoFiles/chem/smiles.csv` open, navigate to **Chem > Calculate >
   Descriptors...** in the top menu. The Descriptors dialog opens.
3. Select at least one descriptor and click **OK**.
4. Observe the outcome within 60 seconds: a visible error notification appears
   (error balloon such as "Timeout of waiting for container status change", or a
   dialog-level error message), the UI does not hang on an indefinite spinner, and
   the table column count is unchanged (no descriptor columns appended).

Expected:

- The error notification is visible within 60 seconds of clicking OK; no indefinite
  spinner or frozen UI.
- The table column count is unchanged — no descriptor columns added on a failed
  container call.
- The error text is human-readable (not a raw stack trace).

---
{
  "order": 5,
  "datasets": ["System:DemoFiles/chem/smiles.csv"]
}
