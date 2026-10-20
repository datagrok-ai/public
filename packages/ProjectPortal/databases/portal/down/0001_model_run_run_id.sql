-- Domain schema portal: 0.2.5 -> 0.2.6 (down)
-- ems-migration: from=8b98daf7bb219ba4f891b2b500b6588cf6fbf5adefd3212ae8802d0684871cf2 to=032f7527cc0526e4831471f4855687567aeef62fc2034d74929fba2c4b41c197
-- Written by hand in the form of the schema snapshot tool (grok schema migrate is not in this checkout).
-- Registry-only changes (no SQL; converged by the next registration):
--   * column-meta model_run.func_call: description

-- add-column model_run.run_id
DROP INDEX IF EXISTS portal.ux_model_run_run_id;
-- Rows recorded by run id alone get their soft reference back from it, so NOT NULL can return.
UPDATE portal.model_run SET func_call = run_id::uuid WHERE func_call IS NULL AND run_id IS NOT NULL;
ALTER TABLE portal.model_run DROP COLUMN IF EXISTS run_id; -- data loss: run ids (kept in func_call above)

-- change-required model_run.func_call: required false -> true [manual]
ALTER TABLE portal.model_run ALTER COLUMN func_call SET NOT NULL;
