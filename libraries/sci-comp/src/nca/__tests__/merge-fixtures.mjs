/**
 * merge-fixtures.mjs — inject the PKNCA moment/lag reference fields produced by
 * regen-fixtures.R (`fixtures/_new_params.json`) into the committed 01/02/03
 * fixtures, preserving every existing value exactly.
 *
 * The existing 8-parameter values are NOT recomputed here — regen-fixtures.R
 * only validates that its PKNCA configuration reproduces them (to ~1e-11). This
 * step purely appends `aumclast, aumcinf_obs, mrt, vss, tlag, pct_aumcextrap`
 * to each profile's `parameters`. Run after regen-fixtures.R; see REGEN.md.
 */
import {readFileSync, writeFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {dirname, join} from 'node:path';

const dir = dirname(fileURLToPath(import.meta.url));
const fixturesDir = join(dir, 'fixtures');
const NEW_KEYS = [
  'aumclast', 'aumcinf_obs', 'mrt', 'vss', 'tlag', 'pct_aumcextrap',
];
/**
 * Fields merged into `provenance` rather than `parameters`. `span_ratio` is a
 * lambda_z FIT DIAGNOSTIC, so it belongs beside the existing `lambda_z_*` keys —
 * not with the reported PK parameters.
 */
const NEW_PROVENANCE_KEYS = ['span_ratio'];

const newParams = JSON.parse(
  readFileSync(join(fixturesDir, '_new_params.json'), 'utf-8'));

for (const dataset of Object.keys(newParams)) {
  const path = join(fixturesDir, `${dataset}.json`);
  const fx = JSON.parse(readFileSync(path, 'utf-8'));
  for (const profile of fx.profiles) {
    const subj = profile.profile_key.subject;
    const np = newParams[dataset][subj];
    if (!np) throw new Error(`${dataset}: no new params for subject ${subj}`);
    for (const k of NEW_KEYS) profile.parameters[k] = np[k];
    if (NEW_PROVENANCE_KEYS.length > 0) {
      profile.provenance = profile.provenance ?? {};
      for (const k of NEW_PROVENANCE_KEYS) profile.provenance[k] = np[k];
    }
  }
  writeFileSync(path, JSON.stringify(fx, null, 2) + '\n', 'utf-8');
  console.log(`merged ${NEW_KEYS.length} parameter + ` +
    `${NEW_PROVENANCE_KEYS.length} provenance fields into ${dataset}.json ` +
    `(${fx.profiles.length} profiles)`);
}
