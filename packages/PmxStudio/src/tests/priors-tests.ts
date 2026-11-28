import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package';
import {PRIORS_FILE, enginePriors} from '../priors';

category('Priors', () => {
  test('every line has a prior', async () => {
    const p = await enginePriors();
    for (const line of ['MC38', 'CT26', 'B16F10'])
      expect(p.lines[line].papers > 20, true, `${line}: ${p.lines[line].papers} papers`);
  });

  // per-paper values stay local (D444): the distributed file holds aggregates per line only
  test('the priors file has no per-paper table', async () => {
    const text = await _package.files.readAsText(PRIORS_FILE);
    expect(/"studies"|"studyid"/.test(text), false, 'a per-paper entry in the distributed priors');
  });

  // Crux stats.tumor_monitor refuses a measurement error of exactly zero (me_c0 = me_c1 = sd_process = 0), and the
  // model's come from this file alone: no line may ship one
  test('every line has a measurement error', async () => {
    for (const [line, l] of Object.entries((await enginePriors()).lines)) {
      const [c0, c1, process] = [l.noise?.c0 ?? 0, l.noise?.c1 ?? 0, l.twin?.sd_process ?? 0];
      expect(c0 >= 0 && c1 >= 0 && process >= 0 && c0 + c1 + process > 0, true,
        `${line}: measurement error ${c0}, ${c1}, process ${process}`);
    }
  });
});
