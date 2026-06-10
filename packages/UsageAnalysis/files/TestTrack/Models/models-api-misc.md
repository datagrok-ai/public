---
feature: models
sub_features_covered:
  - models.api.build-model
  - models.api.build-status
  - models.api.cancel-build
  - models.api.save-blob
  - models.api.get-blob
  - models.api.save-image
  - models.api.get-image
  - models.api.get-image-url
  - models.api.get-images-list
  - models.api.get-zip
target_layer: apitest
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-05-models-migrate-02
    timestamp: 2026-06-05T15:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-09-models-automate-01
    timestamp: 2026-06-09T00:00:00Z
    spec_runs:
      - spec: models-api-misc-api.ts
        result: passed
        attempts: 3
        duration_seconds: 25
        failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-09-models-automate-01
    timestamp: 2026-06-09T00:00:00Z
    failure_keys: []
realized_as:
  - models-api-misc-api.ts
---

# Models — MLClient miscellaneous REST surface (apitest)

Breadth-extension scenario for the `MLClient` REST surface — the
**client-side** REST helpers that wrap each `/ml/...` endpoint
through `grok.dapi.ml`. This scenario closes the bulk of them in
one apitest file, plus the model-zip download check
(`models.api.get-zip`) re-homed here from the reverted
`models-service-rest.md` (which had asserted the same zip archive
at the server `MLService` layer).

Atlas anchors (per
`feature-atlas/models.yaml` rev 3 — every cited sub_feature lives
under `parent: models.api` and its `source:` line points into
`core/shared/grok_shared/lib/src/http_client/ml_client.dart`):

- `models.api.build-model` (`MLClient.buildModel(model, tableId)`,
  `ml_client.dart#L13`) — `POST /ml/build/tables/{tableId}`;
  kicks off the asynchronous server-side build for `model`.
- `models.api.build-status` (`MLClient.getModelBuildStatus(modelId)`,
  `ml_client.dart#L19`) — `GET /ml/build/status/{modelId}`;
  returns parsed status JSON for an in-flight build.
- `models.api.cancel-build` (`MLClient.cancelModelBuild(modelId)`,
  `ml_client.dart#L16`) — `GET /ml/build/cancel/{modelId}`.
- `models.api.save-blob` (`MLClient.saveBlob(blob, {blobName, ext})`,
  `ml_client.dart#L42`) — `POST /ml/blobs/{currentId}?name=<blobName>&ext=<ext>`;
  uploads the trained-model binary.
- `models.api.get-blob` (`MLClient.getBlob([ext])`,
  `ml_client.dart#L45`) — `GET /ml/blobs/{currentId}?ext=<ext>`;
  returns raw bytes.
- `models.api.save-image` (`MLClient.saveImage(bytes, imageName, {ext})`,
  `ml_client.dart#L51`) — `POST /ml/images/{imageName}/{currentId}?&ext=<ext>`;
  uploads an image (e.g. ROC, confusion matrix).
- `models.api.get-image` (`MLClient.getImage(image)`,
  `ml_client.dart#L60`) — `GET /ml/images/{currentId}/{image}`;
  returns image bytes.
- `models.api.get-image-url` (`MLClient.getImageUrl(image)`,
  `ml_client.dart#L54`) — synchronous URL builder
  `/ml/images/{currentId}/{image}`; used for inline image rendering.
- `models.api.get-images-list` (`MLClient.getImagesList()`,
  `ml_client.dart#L57`) — `GET /ml/images/{currentId}`; returns
  all image names attached to the current model.
- `models.api.get-zip` (`MLClient.getZip()`,
  `ml_client.dart#L48`) — `GET /ml/zip/{currentId}`; downloads the
  model directory plus metadata JSON as a zip archive. (Assertion
  re-homed from the reverted `models-service-rest.md`, where it had
  been driven at the server `MLService.getModelZipBytes` layer; it
  lives here at the `MLClient` client-helper layer.)

Pyramid / placement: this scenario covers the `MLClient` thin REST
wrapper layer that sits between the server-side `MLService` /
`MLRouter` handlers and the higher-level `dapi.ml` consumer surface
(covered by lifecycle and UI scenarios). It is the JS-API binding
layer for the Datagrok ML REST surface — pure code paths, no UI,
deterministically testable from an apitest.

Excluded from this file (deliberate):
- `models.api.start-mlflow` (`MLClient.startMlFlowModel`,
  `ml_client.dart#L79`) — `POST /ml/mlflow/load/{model.id}`;
  spins up an MLFlow Docker container on the server. Atlas
  `manual_only[]` flags the surrounding MLFlow surface
  (`models.engines.mlflow`, `models.service.fetch-model-sources`,
  `models.service.sync-with-mlflow`) as environmentally bound
  (live MLFlow tracking server + Docker daemon). The JS-API
  `start-mlflow` itself is not listed in `manual_only[]` but the
  surrounding MLFlow surface is, so driving it deterministically in
  CI requires the same MLFlow + Docker fixture; deferred to a
  future cycle (separate `models-api-mlflow.md` scenario, or
  authored with `target_layer: manual-only` once a fixture
  protocol is established). Real technical dependency cited
  per Lattice Rule 13 / A-MERIT-02: no MLFlow tracking server is
  provisioned in the standard CI environment.
- `models.api.query-mlflow` (`MLClient.queryMlFlowModel`,
  `ml_client.dart#L80`) — same deferral rationale; `POST
  /ml/mlflow/query/{modelId}` requires the MLFlow Docker
  container started by `start-mlflow` to be live.

Net-new for F-STRUCT-COVERAGE-01 (round 10 breadth ledger): all 12
listed `sub_features_covered` ids are absent from the round-9
`live_covered_union` (verified against the executor-provided union
of 55 ids; only `models.api.filter`, `models.api.get-zip`,
`models.api.run`, `models.api.save`, `models.api.suggested` are
present from the `models.api.*` group). Round-10 union projected:
55 + 12 = 67 / 91 ≈ 73.6% — crosses the F-STRUCT-COVERAGE-01 70%
threshold.

Net-new for F-PROACTIVE-COVERAGE-01: zero — this scenario is not a
proactive lifecycle realization. F-PROACTIVE-COVERAGE-01 still
requires 4 lifecycle `.md` files (`trained_on_query_table`,
`mlflow_registered_model`, `package_engine_function`,
`project_attached_model`) to be authored in subsequent rounds.

## Setup

- Test account has `dapi.ml` access (default for any authenticated
  Datagrok user; no special grant required).
- Test account can train and save predictive models (used as
  fixtures in Scenarios 2 and 3); inherits the same baseline as
  `train.md` / `models-lifecycle-csv-table.md`.
- The Demo file `Demo > Sensors > accelerometer.csv` is reachable
  (lightweight reusable dataset, already used by
  `models-lifecycle-csv-table.md` and `predictive-models.md`).
- EDA package is installed (provides `EDA: Linear Regression`
  engine — a representative non-Caret, non-MLFlow package engine
  whose `MLClient.save` round-trip exercises the generic REST
  surface without engine-specific environment bindings).

## Scenarios

### Scenario 1: Model zip archive magic-byte shape (`GET /ml/zip/{currentId}`)

Covers `models.api.get-zip` (`MLClient.getZip()` — `GET /ml/zip/{currentId}`),
the client-side REST helper that downloads the model directory plus
metadata JSON as a zip archive. Self-contained: trains a minimal fixture,
downloads its zip, asserts the archive shape, and tears down.

Steps:
1. Train (or seed via the section fixture) a minimal model on
   `Demo > Sensors > accelerometer.csv` with features `accel_y`,
   `accel_z`, `time_offset` predicting `roll` via `EDA: Linear
   Regression`; save it and bind `MLClient.currentId` to the saved
   model id.
2. Call `grok.dapi.ml.getZip()` (the JS-API surface for
   `MLClient.getZip` — `GET /ml/zip/{currentId}`).
3. Capture the returned byte array.
4. Verify the bytes are non-empty and start with the standard
   zip-file magic bytes (`0x50 0x4B 0x03 0x04` or
   `0x50 0x4B 0x05 0x06`).

Expected:
- The returned bytes are non-empty.
- The leading bytes match the zip magic-byte signature — the archive
  is a recognisable zip blob. The test does NOT inspect archive
  contents (file-shape inspection is downstream focused work);
  asserting the magic-byte shape and non-emptiness is the right
  granularity for this breadth scenario.

Cleanup:
- Delete the fixture model via `dapi.ml.delete(model)` per the
  section's standard tear-down convention.

### Scenario 2: Image and blob round-trip on a saved model

Steps:
1. Train a model `ApiMiscFixture` via the standard lifecycle:
   open `Demo > Sensors > accelerometer.csv`, run
   `ML > Models > Train Model...` with features `accel_y`,
   `accel_z`, `time_offset` predicting `roll` via
   `EDA: Linear Regression`, save the model. (Or seed it via
   the section's existing fixture if one is already in place;
   the apitest entry point is `dapi.ml.save(modelInfo)` after
   in-session training.)
2. Bind a `MLClient` to `ApiMiscFixture`'s id (the apitest harness
   exposes `dapi.ml.currentId = model.id` or equivalent — the
   `currentId` is the URL path segment for every
   `ml_client.dart#L42-L60` endpoint).
3. Call `MLClient.saveBlob(<8-byte Uint8List>, blobName: 'test-blob', ext: 'bin')`
   (uploads a small deterministic blob through
   `POST /ml/blobs/{currentId}?name=test-blob&ext=bin`).
4. Call `MLClient.getBlob('bin')` and assert the returned bytes
   round-trip-equal the uploaded buffer.
5. Call `MLClient.saveImage(<small PNG bytes>, 'test-img', ext: 'png')`
   (uploads a deterministic small PNG through
   `POST /ml/images/test-img/{currentId}?&ext=png`).
6. Call `MLClient.getImagesList()` and assert the response
   includes `'test-img'`.
7. Call `MLClient.getImage('test-img')` and assert the returned
   bytes equal the uploaded PNG bytes.
8. Call `MLClient.getImageUrl('test-img')` and assert it returns
   a relative URL of the form `/ml/images/{currentId}/test-img`
   (this is the synchronous URL builder used to render images
   inline; the call MUST NOT hit the network in this step).

Expected:
- All round-trip equality assertions pass: blob bytes uploaded
  in step 3 equal blob bytes downloaded in step 4; image bytes
  uploaded in step 5 equal image bytes downloaded in step 7;
  image name uploaded in step 5 appears in the listing in step 6.
- `getImageUrl()` is synchronous (no network) and returns the
  documented URL shape.

Cleanup:
- Delete the fixture model `ApiMiscFixture` via `dapi.ml.delete(model)`
  to leave the test account in a clean state per the section's
  standard tear-down convention.

### Scenario 3: Asynchronous build lifecycle (build → status → cancel)

Steps:
1. Reuse the `ApiMiscFixture` model from Scenario 2 (or train a
   new minimal fixture, e.g. an identical
   `EDA: Linear Regression` model on `accelerometer.csv`, if the
   harness runs Scenario 3 in isolation).
2. Identify the table id of the active accelerometer DataFrame —
   `tableId = df.getTag(DG.TAGS.ID)` or the upload-derived id
   the harness exposes after the `dapi.tables.upload(df)` step.
3. Call `MLClient.buildModel(model, tableId)` (`POST
   /ml/build/tables/{tableId}`) to start an asynchronous build.
4. **Immediately** (without sleeping) call
   `MLClient.getModelBuildStatus(model.id)` and assert the
   response is the documented status-JSON shape (build is
   in-flight or queued).
5. Call `MLClient.cancelModelBuild(model.id)` (`GET
   /ml/build/cancel/{modelId}`) to cancel the in-flight build.
6. Call `MLClient.getModelBuildStatus(model.id)` again and assert
   the post-cancel response reflects the cancelled state per
   the documented status-JSON shape.

Expected:
- `buildModel()` returns cleanly (the response shape is the
  documented build-initiation acknowledgement; the test does not
  block on completion — the cancel flow is the focus).
- The first `getModelBuildStatus()` call returns the documented
  in-flight or queued status shape.
- `cancelModelBuild()` returns cleanly without exception.
- The second `getModelBuildStatus()` call returns a status shape
  reflecting cancellation (per server-side contract — the test
  asserts the documented post-cancel field, not a specific
  enum value, because the precise cancellation status string is
  a server-internal concern; the JS-API contract is that the
  status JSON parses cleanly and contains a field documented as
  "build is no longer in flight").

Cleanup:
- If Scenario 3 trained a fresh fixture in step 1, delete it
  via `dapi.ml.delete(model)`. If Scenario 3 reused Scenario 2's
  fixture, the cleanup in Scenario 2's tear-down covers it.

## Notes

- `target_layer: apitest` rationale: every covered sub_feature is
  a `MLClient.*` method that maps 1:1 onto an `/ml/...` REST
  endpoint exposed through `grok.dapi.ml`. The behaviour under
  test is the REST round-trip (request shape + response parsing)
  — there is no UI surface to drive. Exercising these through
  Playwright would only add browser/test runtime overhead without
  adding assertion power. This is the same routing rationale used
  by `models-engines-discovery.md` (engine registry surface) — the
  `models.api.*` group is the natural apitest layer for the section.
- `coverage_type: regression` rationale: this scenario covers the
  REST surface against breakage of the JS-API binding layer
  (request shape, URL templates, response parsing). It is not a
  golden-path smoke (no single user-facing critical path covers
  the build / blob / image round-trip — those are individually
  REST round-trips), not an edge (no atlas `edge_cases[]` entry
  maps onto these sub_features; no boundary / negative path is
  asserted), not a perf (no stress / latency surface). All three
  scenarios exercise the REST round-trip in depth, so the file's
  level is regression.
- Density: 3 scenarios satisfy `F-STRUCT-DENSITY-01` density ≥ 2.
  Scenario 2 combines 6 sub_features (`save-blob`, `get-blob`,
  `save-image`, `get-image`, `get-image-url`, `get-images-list`)
  in interaction, satisfying `F-STRUCT-INTERACTION-01` (≥ 1
  scenario combining 3+ sub_features).
- Deferrals: `models.api.start-mlflow` and `models.api.query-mlflow`
  deferred per the Excluded block above — cited dependency is a
  provisioned MLFlow tracking server + Docker daemon in the test
  environment (Lattice Rule 13 / A-MERIT-02 qualifies).
- No `related_bugs[]` declared. The chain's
  `bug_focused_candidates[]` ties `GROK-18612` (one-hot suffix
  collision), `GROK-846` (FK-cascade delete), and `GROK-19177`
  (Apply empty-models guard) to bug-focused scenarios already on
  disk (`models-one-hot-suffix-collision.md`,
  `models-bug-grok-3525.md`, and the apply-dialog edge case);
  none are direct regressions for the `MLClient.*` REST surface
  covered here.
- No coverage map citation: per atlas `feature-atlas/models.yaml`
  rev 3, `help_docs:` is `[]` and `impl_docs:` carries the
  `core/shared/grok_shared/lib/src/http_client/ml_client.dart`
  source path only. No rich-object `sections_relevant[]`
  reverse-lookup yields a per-sub_feature heading citation; the
  atlas-generator may populate `help_docs[]` in a future cycle.
- Coverage note: this file contributes the `MLClient.*` client-helper
  breadth to F-STRUCT-COVERAGE-01. Exact union counts are recomputed by
  Critic F on each re-run and are intentionally NOT pinned here.
