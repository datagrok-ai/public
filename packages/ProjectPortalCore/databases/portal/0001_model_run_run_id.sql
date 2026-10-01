-- Domain schema portal: 0.2.5 -> 0.2.6
-- ems-migration: from=032f7527cc0526e4831471f4855687567aeef62fc2034d74929fba2c4b41c197 to=8b98daf7bb219ba4f891b2b500b6588cf6fbf5adefd3212ae8802d0684871cf2
-- Written by hand in the form of the schema snapshot tool (grok schema migrate is not in this checkout).
-- Registry-only changes (no SQL; converged by the next registration):
--   * column-meta model_run.func_call: description

-- change-required model_run.func_call: required true -> false [manual]
-- The server cannot write the soft reference to Core.func_calls yet (soft-ref alias bug), so a run
-- is recorded by run_id alone until it can. Migration debt: backfill func_call from run_id then.
ALTER TABLE portal.model_run ALTER COLUMN func_call DROP NOT NULL;

-- add-column model_run.run_id
ALTER TABLE portal.model_run ADD COLUMN IF NOT EXISTS run_id text;
CREATE UNIQUE INDEX IF NOT EXISTS ux_model_run_run_id ON portal.model_run (run_id) WHERE NOT is_deleted;
