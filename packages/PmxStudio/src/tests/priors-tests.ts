import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {enginePriors, mapPrior} from '../priors';

// The prior file was computed offline by Crux's native CLI;
// the package recomputes it with Crux's wasm module from the same per-paper table: the two must agree
// (stats.map_prior is deterministic quadrature, declared exact to 1e-6).
category('Priors', () => {
  test('map_prior wasm = native', async () => {
    const p = await enginePriors();
    for (const [line, l] of Object.entries(p.lines)) {
      const r = await mapPrior(line);
      const mean = r.main.col['map_mean'][0];
      const sd = r.main.col['map_sd'][0];
      expect(Math.abs(mean - l.log_alpha.map_mean) < 1e-6, true, `${line} map_mean ${mean} vs ${l.log_alpha.map_mean}`);
      expect(Math.abs(sd - l.log_alpha.map_sd) < 1e-6, true, `${line} map_sd ${sd} vs ${l.log_alpha.map_sd}`);
    }
  }, {timeout: 120000});

  test('every line has a prior', async () => {
    const p = await enginePriors();
    for (const line of ['MC38', 'CT26', 'B16F10'])
      expect(p.lines[line].studies.length > 20, true, `${line}: ${p.lines[line].studies.length} studies`);
  });
});
