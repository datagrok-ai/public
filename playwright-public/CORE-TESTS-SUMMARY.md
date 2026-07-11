# Core playwright-public — CI Green-Up Summary

**Branch:** `opavlenko/playwright-public-ci-fixes-0709` (not merged to master) · **Date:** 2026-07-10
**Scope:** the **core / platform** suites that stay in `public/playwright-public/` (not the package-owned suites — those are in `PACKAGE-PLAYWRIGHT-CODE-FINDINGS.md`).
**Job:** Jenkins **Test-Playwright** (`PLAYWRIGHT_DIR=public/playwright-public`, per-suite `TEST_GREP`). Read the real result from the `.jest` artifact (pipeline SUCCESS ≠ tests passed).

## Result: all 14 core suites GREEN (0 failures)

Time = pure Playwright run time (per-suite `.jest`, excludes prereq package builds).

| Suite | Passed | Failed | Skipped | Time | Status | Fix |
|---|---:|---:|---:|---:|---|---|
| scripts | 29 | 0 | 2 | 4:03 | ✅ Green | — |
| spaces | 18 | 0 | 0 | 4:23 | ✅ Green | — |
| user groups | 20 | 0 | 0 | 6:02 | ✅ Green | — |
| stickymeta | 3 | 0 | 0 | 1:45 | ✅ Green | — |
| projects | 38 | 0 | 1 | 12:15 | ✅ Green | deep-link reopen via JS API (`project.open()`, not cold `/p/` goto) |
| Sharing | 7 | 0 | 0 | 7:27 | ✅ Green | — |
| Chat | 1 | 0 | 0 | 2:26 | ✅ Green | author label tolerant to login vs friendlyName across sessions |
| Models | 7 | 0 | 0 | 3:58 | ✅ Green | — |
| Viewers | 27 | 0 | 0 | 16:38 | ✅ Green | pivot re-select; trellis menu via dispatch + DevTools-conditional; Map/annotation/forms/grid via prereq |
| Viewers/FilterPanel | 4 | 0 | 0 | 1:28 | ✅ Green | — |
| legend | 16 | 0 | 0 | 10:35 | ✅ Green | `tv.loadLayout()` after `project.open()` (reopen restores Grid but not custom viewers) ×4 |
| browse | 98 | 0 | 15 | 18:57 | ✅ Green | Filter-03: assert a real Apps quick-filter ("Favorites"), not the non-existent "Used by me" |
| connections | 21 | 0 | 12 | 4:22 | ✅ Green | correct datagrok Postgres pw in Jenkins vault + CasC reload; 09 write-queries skipped |
| queries | 12 | 0 | 7 | 3:36 | ✅ Green | same datagrok Postgres pw (Jenkins reload) |
| **TOTAL** | **301** | **0** | **37** | **~1h55m** | **✅ Green** | |

## What was fixed (test-harness)

- **projects** — cold `page.goto('/p/<ns>.<name>')` did not open a just-saved project (SPA showed Home; the by-path route lags the search index). Reopen via `dapi.projects.find(id).open()` and assert id + rowCount + URL.
- **legend ×4 / Charts-radar (package)** — `project.open()` re-materializes the table + default Grid but **not** the attached custom-viewer layout; load the saved layout explicitly after reopen.
- **Chat** — the comment author label renders as the login (`test.user`) in one session and the friendlyName (`Test`) in another; assert against the identity set {login, friendlyName, captured} rather than one exact string.
- **Viewers** — pivot-table re-select after rowSource switch; trellis "To Script" driven via a dispatched `contextmenu` and made conditional on the DevTools plugin (it contributes the menu item and can't build on the CI stack); Map/annotation/forms/grid fixed by adding the GIS/PowerPack/PowerGrid prereq packages.
- **browse Filter-03** — the Apps gallery has only All / Favorites / Created recently (no "Used by me" — that belongs to the connections/dockerfiles galleries), so assert "Favorites".

## Infrastructure

- **connections/queries creds** — the 4 failures were a wrong `datagrok` Postgres password in the Jenkins credential store. Rotated the age-encrypted `test_postgres_datagrok_user__password` in `infra/ansible/jenkins-credentials/secrets.yaml`; the pipeline binds it as `DG_PG_*`. Fixed only after a **Jenkins CasC restart** (credentials load at startup — a config redeploy alone doesn't reload them). Verified the password directly against `db.datagrok.ai:54322`.
- **connections/09** external-provider write-queries (CREATE/INSERT/UPDATE/DROP) **skipped** on CI — needs a writable external Postgres; the CI stack has none (datagrok is read-only on Northwind; the `:54327` test-DB credentials are unavailable). Agreed.
- **Jenkins job access** — added a CasC role-strategy item role granting `opavlenko` Job/Read+Build+Cancel on the `dg-fix-reports` job (it had been set in the UI and got wiped on a CasC reload; codified in `jenkins.yaml`).
- **Pipeline note** — Test-Playwright runs on the `stress` worker node; after a Jenkins restart that node must be brought back online or builds hang with "There are no nodes with the label 'stress'".

## Notes

- **37 skips** are all pre-existing static `test.skip` (e.g. scripts Pyodide/standalone, projects cross-namespace move), absent packages/backends on the minimal CI stack (browse ×15), or the agreed connections/09 write-queries. Nothing newly skipped to hide a failure.
- Per-suite source builds (Test-Playwright): scripts/spaces/user-groups #230; stickymeta/Models #236; projects/Sharing #242; Chat #246; Viewers #241+#249; FilterPanel/legend #241; browse #237+#240+#247; connections/queries #251.
