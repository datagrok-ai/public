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
  it was never provisioned). Then open System:DemoFiles/chem/smiles.csv and go to
  Chem > Calculate > Descriptors... The failure surfaces at dialog OPEN, not at OK:
  the descriptor tree is fetched from the container while the dialog builds, so with
  the container down the dialog opens with the message "Could not load descriptors.
  The Chem service may be unavailable." where the tree should be, and the OK button
  never enables (no tree means nothing is selected, so the dialog stays invalid —
  descriptors-calculation.ts, initTree/isValid). Confirm that message appears within
  60 seconds instead of an endless spinner, that OK cannot be pressed, and that no
  descriptor columns were appended (column count unchanged).
automation_candidate: false
blocked_by: []
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

**Where the failure surfaces: at dialog open, not at OK.** The descriptor tree is fetched
from the `chem-chem` container while the dialog is being built. With the container down
that fetch fails, the dialog shows a message where the tree should be, and OK never
becomes clickable — an empty tree means nothing is selected, so the dialog never becomes
valid. Do not wait for a balloon after pressing OK: you will never get to press OK, and
recording "no balloon after OK" would report a failure that is really this earlier,
correct symptom.

Steps:

1. Confirm the `chem-chem` container is not running.
2. With `System:DemoFiles/chem/smiles.csv` open, note the table's current column count
   (the grid's column headers, or the table's tooltip in **Tables**). Write it down —
   the last step compares against it.
3. Navigate to **Chem > Calculate > Descriptors...** in the top menu. The Descriptors
   dialog still opens: the **Table** and **Molecules** inputs are present and populated.
4. Look at the area where the descriptor tree normally is, and start counting from the
   moment the dialog opened. Within 60 seconds it must show the text **"Could not load
   descriptors. The Chem service may be unavailable."** — no tree, no endless spinner.
5. Try the **OK** button. It must be disabled and must not submit anything. Cancel the
   dialog.
6. Re-check the table's column count against the number from step 2. It must be
   unchanged.

Expected:

- The dialog resolves to a readable message within 60 seconds of opening; no indefinite
  spinner and no frozen UI. This bounded resolution is the GROK-17621 guard.
- OK stays disabled while the tree is unavailable, so no container call is ever submitted.
- The table column count is unchanged — no descriptor columns are added.

What to record, case by case:

- **Message appears within 60 s, OK disabled, column count unchanged.** PASS. Record the
  message text verbatim and roughly how long it took to appear.
- **Message never appears — spinner or blank tree area past 60 s.** FAIL, and this is the
  regression GROK-17621 guards. Record how long you actually waited before giving up.
- **A timeout balloon instead (e.g. "Timeout of waiting for container status change").**
  Still a PASS on boundedness if it arrives within 60 s and OK never submitted a call, but
  it is a different symptom than the one above — record which one you saw, because it
  means the container was registered and starting rather than absent.
- **OK is enabled, or the descriptor tree loads normally.** NOT a result either way: the
  container was reachable, so the precondition was never established. Stop, take the
  container down properly, and start again. Do not record a pass or a fail.
- **A raw stack trace or minified-JS text instead of the message.** FAIL. The error
  surfacing is only half the contract; an unreadable error still fails the step. (A stack
  trace logged to the browser console is expected and is not this failure — the failing
  fetch is logged there deliberately. Judge only what is shown in the dialog.)
- **Descriptor columns appear in the table.** The container served the request; treat it
  as a precondition failure, as above, not as a pass.

## Cleanup

Nothing is created by this scenario, so there is nothing to delete. Restore the
environment you changed: if you stopped the `chem-chem` container to establish the
precondition, start it again, and confirm it reaches a running state before leaving
the stand for anyone else.

---
{
  "order": 5,
  "datasets": ["System:DemoFiles/chem/smiles.csv"]
}
