import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {_package} from '../package-test';
import {readDataframe} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, after, expect, category, test, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {cliffMcs, createPropPanelElement, createTooltipElement} from '../analysis/activity-cliffs';
import {MALFORMED_DATA_WARNING_CLASS} from '../constants';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import { PackageFunctions } from '../package';
import {ChemicalEngine, cruxReady, setChemicalEngine} from '../engine/chem-engine';
import {cruxFallbackCount} from '../engine/fallback';
import {getCruxService} from '../crux/crux-searches';
import {getMCS} from '../utils/most-common-subs';
import {getUncommonAtomsAndBonds} from '../utils/chem-common';

/**
 * Cliff pairs of the tests' data (file, column, rows) and what RDKit's Rascal MCES leaves out of each molecule (atoms,
 * bonds), as Crux mode highlights it: the pinned RDKit's `rdRascalMCES.FindMCES` at Chem's options (similarityThreshold
 * 0.3, maxBondMatchPairs 10000, timeout 5). The triterpenoid pair needs maxBondMatchPairs over RDKit's 1000. The last
 * pair's MCES is not its FMCS: 12 bonds in two fragments (Rascal leaves out partly matched aromatic rings) where FMCS
 * finds 18 connected ones.
 */
const RASCAL_PAIRS: [string, string, number, number, number[][][]][] = [
  ['tests/activity_cliffs_test.csv', 'smiles', 4, 5, [[[0], [0]], [[12], [11]]]],
  ['tests/activity_cliffs_test.csv', 'smiles', 10, 11, [[[], []], [[], []]]],
  ['tests/spgi-100.csv', 'Structure', 32, 42, [[[0, 1], [0, 1]], [[0, 32], [0, 31]]]],
  ['tests/Test_smiles_malformed.csv', 'canonical_smiles', 8, 15, [[[], []], [[17], [16]]]],
  ['tests/activity_cliffs_test.csv', 'smiles', 4, 9, [
    [[0, 1, 2, 3, 4, 5, 6, 17, 18, 19, 20, 21, 24, 25],
      [0, 1, 2, 3, 4, 5, 6, 8, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26]],
    [[7, 8, 9, 10, 11, 12, 13, 14, 15], [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 21, 23, 24]]]],
];

/** A cliff of two 160-atom peptides (6,904 matching bond pairs; their FMCS did not finish in 13 minutes): Crux mode
 * only. */
const PEPTIDE_PAIR: [string, string, number, number, number[][][]] = ['v3000_sample.csv', 'molecule', 0, 4, [
  [[52, 120, 121], [52, 56, 116, 122, 123]], [[52, 53, 119, 120, 121], [52, 53, 57, 117, 121, 122, 123, 124]]]];

/** The molecules of a {@link RASCAL_PAIRS} pair. */
async function pairOf(file: string, column: string, i: number, j: number): Promise<[string, string]> {
  const col = (await readDataframe(file)).col(column)!;
  return [col.get(i), col.get(j)];
}

/** The pairs of `molecules` whose Morgan (radius 2, 2048 bits) Tanimoto similarity reaches `cutoff`: its cliffs'. */
function similarPairs(molecules: string[], cutoff: number): [number, number][] {
  const rdkit = chemCommonRdKit.getRdKitModule();
  const fps = molecules.map((m) => {
    const mol = m ? rdkit.get_mol(m) : null;
    try {
      return mol ? mol.get_morgan_fp(JSON.stringify({radius: 2, nBits: 2048})) : null;
    } finally {
      mol?.delete();
    }
  });
  const pairs: [number, number][] = [];
  for (let i = 0; i < fps.length; i++) {
    for (let j = i + 1; j < fps.length; j++) {
      const [a, b] = [fps[i], fps[j]];
      if (!a || !b)
        continue;
      let both = 0;
      let either = 0;
      for (let k = 0; k < a.length; k++) {
        both += +(a[k] === '1' && b[k] === '1');
        either += +(a[k] === '1' || b[k] === '1');
      }
      if (either && both / either >= cutoff)
        pairs.push([i, j]);
    }
  }
  return pairs;
}

const rdkitPoolAsks = () => chemCommonRdKit._rdKitService?.timesInitialized ?? 0;

/** Median, 95th percentile and maximum of `ms`. */
function timings(ms: number[]): string {
  const sorted = [...ms].sort((a, b) => a - b);
  const at = (q: number) => sorted[Math.min(sorted.length - 1, Math.floor(q * sorted.length))].toFixed(1);
  return `median ${at(0.5)}, p95 ${at(0.95)}, max ${sorted[sorted.length - 1].toFixed(1)} ms`;
}


category('top menu activity cliffs', async () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('activityCliffsOpenAndLayoutApply.smiles', async () => {
    const df = DG.Test.isInBenchmark ?
      await readDataframe('tests/smiles_1K_with_activities.csv') :
      await readDataframe('tests/activity_cliffs_test.csv');
    await _testActivityCliffsOpen(df, 'smiles', 'Activity', DG.Test.isInBenchmark ? 60 : 80,
      DG.Test.isInBenchmark ? 16 : 2, !DG.Test.isInBenchmark);
  }, {timeout: 20000, benchmark: true, benchmarkTimeout: 180000});

  test('activityCliffs.cutoff20', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/activity_cliffs_test.csv'), 'smiles', 'Activity', 20, 52);
  });

  test('activityCliffs.cutoff95', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/activity_cliffs_test.csv'), 'smiles', 'Activity', 95, 1);
  });

  test('activityCliffsOpen.molV2000', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/spgi-100.csv'), 'Structure', 'Chemical Space X', 80, 1);
  }, {timeout: 20000});

  test('activityCliffsOpen.molV3000', async () => {
    await _testActivityCliffsOpen(await readDataframe('v3000_sample.csv'), 'molecule', 'Activity', 80, 185);
  });

  test('activityCliffs.emptyValues', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/activity_cliffs_empty_rows.csv'),
      'smiles', 'Activity', 80, 1);
  });

  test('activityCliffs.malformedData', async () => {
    DG.Balloon.closeAll();
    await _testActivityCliffsOpen(await readDataframe('tests/Test_smiles_malformed.csv'),
      'canonical_smiles', 'FractionCSP3', 80, 24);
    try {
      await awaitCheck(() => {
        console.log(document.querySelector(`.${MALFORMED_DATA_WARNING_CLASS}`)?.innerHTML);
        return document.querySelector(`.${MALFORMED_DATA_WARNING_CLASS}`)?.innerHTML ===
        '3 molecules with indexes 14,31,41 are possibly malformed and are not included in analysis';
      }, 'cannot find warning balloon', 5000);
    } finally {
      grok.shell.closeAll();
      DG.Balloon.closeAll();
    }
  });

  test('activityCliffs.viewReopenedDuringRun', async () => {
    const df = await readDataframe('tests/activity_cliffs_test.csv');
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    const encodingFunc = DG.Func.find({name: 'getFingerprints', package: 'Chem'})[0];
    try {
      const run = PackageFunctions.activityCliffs(df, df.col('smiles')!, df.col('Activity')!, 80,
        DimReductionMethods.UMAP, BitArrayMetricsNames.Tanimoto, encodingFunc, undefined, undefined, true);
      tv.close();
      const reopened = df.clone();
      reopened.name = df.name;
      const reopenedView = grok.shell.addTableView(reopened);
      await run;
      expect(grok.shell.tv.dataFrame === df, true);
      expect(Array.from(reopenedView.viewers).some((v) => v.type === DG.VIEWER.SCATTER_PLOT), false);
      expect(Array.from(grok.shell.tv.viewers).some((v) => v.type === DG.VIEWER.SCATTER_PLOT), true);
    } finally {
      grok.shell.closeAll();
    }
  });

  test('activityCliffs_layout', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/spgi-100.csv'), 'Structure', 'Chemical Space X', 80, 1);
  });

  test('cliffMcs.crux', async () => {
    // Crux mode: what the tooltip highlights is what Rascal's MCES (the Crux workers) leaves out; no RDKit worker asked
    setChemicalEngine(ChemicalEngine.Crux);
    try {
      expect(await cruxReady(), true, 'the Crux module loads');
      await getCruxService();
      const [asks, fallbacks] = [rdkitPoolAsks(), cruxFallbackCount('activity cliffs MCS')];
      for (const [file, column, i, j, uncommon] of [...RASCAL_PAIRS, PEPTIDE_PAIR]) {
        const mcs = await cliffMcs(await pairOf(file, column, i, j));
        const got = Array.isArray(mcs) ? mcs.map((s) => [s.atoms, s.bonds]) : mcs;
        expect(JSON.stringify(got), JSON.stringify(uncommon), `${file} rows ${i}, ${j}`);
      }
      expect(cruxFallbackCount('activity cliffs MCS'), fallbacks, 'RDKit found an MCS for Crux');
      expect(rdkitPoolAsks(), asks, 'the RDKit workers were asked for');
    } finally {
      setChemicalEngine(null);
    }
  }, {timeout: 60000});

  test('cliffMcs.rdkit', async () => {
    // RDKit mode: FMCS's SMARTS from the RDKit workers, as before; its highlight is not Rascal's on the last pair
    setChemicalEngine(ChemicalEngine.RDKit);
    try {
      for (const [file, column, i, j] of RASCAL_PAIRS) {
        const pair = await pairOf(file, column, i, j);
        expect(await cliffMcs(pair), await getMCS(DG.Column.fromStrings('smiles', pair), true, true),
          `${file} rows ${i}, ${j}`);
      }
      const [file, column, i, j, uncommon] = RASCAL_PAIRS[RASCAL_PAIRS.length - 1];
      const pair = await pairOf(file, column, i, j);
      const rdkit = chemCommonRdKit.getRdKitModule();
      const fmcs = rdkit.get_qmol(await cliffMcs(pair) as string);
      try {
        const highlighted = pair.map((m) => getUncommonAtomsAndBonds(m, fmcs, rdkit)).map((s) => [s!.atoms, s!.bonds]);
        expect(JSON.stringify(highlighted) !== JSON.stringify(uncommon), true, 'FMCS highlights what Rascal does');
      } finally {
        fmcs.delete();
      }
    } finally {
      setChemicalEngine(null);
    }
  });

  test('cliffMcs.everyPair', async () => {
    // every cliff of the tests' data and of the demo's gets an MCES at Chem's options (with RDKit's default
    // maxBondMatchPairs, 22 of the 24 triterpenoid cliffs got none); and the time per pair, Rascal on a Crux worker and
    // FMCS on the RDKit workers, each a round trip from here (FMCS not on the triterpenoids: 13 to 88 s a pair)
    const times: {[engine: string]: number[]} = {FMCS: [], Rascal: []};
    const missing: string[] = [];
    const service = await getCruxService();
    for (const [file, column, cutoff, fmcs] of [['tests/activity_cliffs_test.csv', 'smiles', 0.2, true],
      ['demo_files/sar_small.csv', 'smiles', 0.78, true], ['tests/Test_smiles_malformed.csv', 'canonical_smiles', 0.8,
        false]] as const) {
      const molecules = (await readDataframe(file)).col(column)!.toList() as string[];
      const pairs = similarPairs(molecules, cutoff);
      console.log(`top menu activity cliffs | cliffMcs.everyPair: ${file}, ${pairs.length} pairs`);
      for (const [i, j] of pairs) {
        let start = performance.now();
        if (!await service.pairUncommonAtomsAndBonds(molecules[i], molecules[j]))
          missing.push(`${file} rows ${i}, ${j}`);
        times.Rascal.push(performance.now() - start);
        if (fmcs) {
          start = performance.now();
          await getMCS(DG.Column.fromStrings('smiles', [molecules[i], molecules[j]]), true, true);
          times.FMCS.push(performance.now() - start);
        }
      }
    }
    const result = `Rascal ${times.Rascal.length} pairs, ${timings(times.Rascal)}; FMCS ${times.FMCS.length} pairs, ` +
      timings(times.FMCS);
    console.log(`top menu activity cliffs | cliffMcs.everyPair: ${result}`);
    expect(missing.join('; '), '', 'cliffs without an MCES');
    return result;
  }, {timeout: 300000}); // the RDKit workers may first finish a long FMCS of an earlier test's cliffs panel
});

async function _testActivityCliffsOpen(df: DG.DataFrame, molCol: string, activityCol: string, simCutOff: number,
  numberCliffs: number, layout?: boolean) {
  await grok.data.detectSemanticTypes(df);
  const actCliffsTableView = grok.shell.addTableView(df);
  await awaitCheck(() => actCliffsTableView.name?.toLowerCase() === df.name?.toLowerCase(),
    'Activity cliffs table view hasn\'t been created', 1000);
  if (molCol === 'molecule') actCliffsTableView.dataFrame.rows.removeAt(51, 489);
  const encodingFunc = DG.Func.find({name: 'getFingerprints', package: 'Chem'})[0];
  await PackageFunctions.activityCliffs(df, df.col(molCol)!, actCliffsTableView.dataFrame.getCol(activityCol), simCutOff,
    DimReductionMethods.UMAP, BitArrayMetricsNames.Tanimoto, encodingFunc, undefined, undefined, true);
  let scatterPlot: DG.Viewer | null = null;
  for (const i of actCliffsTableView.viewers) {
    if (i.type == DG.VIEWER.SCATTER_PLOT)
      scatterPlot = i;
  }
  expect(scatterPlot != null, true);
  const checkScatterPlotInitialized = () => {
    const cliffsLink = scatterPlot!.root.getElementsByClassName('scatter_plot_link');
    return cliffsLink.length > 0 && ((cliffsLink[0] as HTMLElement).innerText.toLowerCase() === `${numberCliffs} cliffs`)
  }
  await awaitCheck(() => {
    return checkScatterPlotInitialized();
  }, 'activity cliffs failed', layout ? 2000 : 60000);
  if (layout) {
    const layout = actCliffsTableView.saveLayout();
    scatterPlot?.close();
    await delay(100);
    actCliffsTableView.loadLayout(layout);
    await delay(100); //to add GridViewer
    //waiting for layout to be applied
    await awaitCheck(() => {
      return checkScatterPlotInitialized();
    }, 'activity cliffs failed', 5000);
  }

  actCliffsTableView.close();
}
