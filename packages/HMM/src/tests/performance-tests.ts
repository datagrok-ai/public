import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {HmmerPool} from '../pool';
import {allowedChains, extractSequence} from '../numbering';
import ANARCI_FIXTURE from './anarci-fixture.json';

/** Generous in-browser budgets (CI machines vary); the tight speed and size
 * budgets are enforced by the Rusty-HMMER harness (web/budgets.json). */
const COLD_START_MS = 5_000;
const MS_PER_CHAIN = 30;

category('HMMER performance', () => {
  test('cold start and numbering throughput', async () => {
    HmmerPool.get(_package).terminate();
    const sequence = ANARCI_FIXTURE[0].sequence;
    let start = performance.now();
    await HmmerPool.get(_package).anarci([['0', sequence]], {scheme: 'imgt', assignGermline: true}, 1);
    const cold = performance.now() - start;

    const df: DG.DataFrame = await grok.dapi.files.readCsv('System:AppData/Bio/samples/antibodies.csv');
    const chains: [string, string][] = [];
    for (const name of ['AntibodyHC', 'AntibodyLC']) {
      const col = df.col(name);
      if (!col) continue;
      for (let i = 0; i < df.rowCount; i++) chains.push([String(chains.length), extractSequence(col.get(i))]);
    }
    start = performance.now();
    const results = await HmmerPool.get(_package).anarci(chains,
      {scheme: 'imgt', allow: allowedChains('imgt'), assignGermline: true});
    const elapsed = performance.now() - start;
    const numbered = results.filter((r) => r.numbered).length;
    const perChain = elapsed / chains.length;
    grok.shell.info(`HMMER: cold start ${cold.toFixed(0)} ms; ${chains.length} chains in ${elapsed.toFixed(0)} ms ` +
      `(${perChain.toFixed(2)} ms/chain, ${HmmerPool.maxWorkers} workers), ${numbered} numbered`);
    console.log(JSON.stringify({coldStartMs: cold, chains: chains.length, elapsedMs: elapsed, perChain,
      workers: HmmerPool.maxWorkers, numbered}));
    expect(cold < COLD_START_MS, true, `cold start ${cold.toFixed(0)} ms`);
    expect(perChain < MS_PER_CHAIN, true, `${perChain.toFixed(2)} ms per chain`);
    expect(numbered > chains.length * 0.9, true, `${numbered} of ${chains.length} numbered`);
  }, {timeout: 300_000});
});
