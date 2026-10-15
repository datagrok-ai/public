import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, delay, expect, expectFloat, test, testViewer}
  from '@datagrok-libraries/test/src/test';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {_package} from '../package-test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {highlightMonomerPosition, scaleActivity} from '../utils/misc';
import {startAnalysis} from '../widgets/peptides';
import {MonomerPosition, MostPotentResidues, SAR_PROPERTIES, SARViewer, SELECTION_MODE} from '../viewers/sar-viewer';
import {COLUMNS_NAMES, SCALING_METHODS, SEM_TYPES} from '../utils/constants';
import {CLUSTER_TYPE, LogoSummaryTable, LST_PROPERTIES} from '../viewers/logo-summary';
import {MutationCliffsViewer} from '../viewers/mutation-cliffs-viewer';
import {PositionHeight} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {analyzeDemoPeptides, DEMO_COLUMN_NAMES, loadDemoPeptides, singleMutationPairs, splitSequences,
  TEST_COLUMN_NAMES} from './utils';
import {showTooltip} from '../utils/tooltips';
import {PeptideUtils} from '../peptideUtils';
import {GEN_COLUMN_NAMES, PEPTIDE_GENERATION_PROPS, PEPTIDE_GENERATION_VIEWER_NAME,
  PeptideGenerationViewer} from '../viewers/peptide-generation-viewer';
import {ACTIVITY_TARGET} from '../utils/constants';

// category('Viewers: Basic', () => {
//   let df: DG.DataFrame;

//   before(async () => {
//     await PeptideUtils.loadComponents();

//     df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
//     const clusterCol = df.col('cluster')!.convertTo(DG.COLUMN_TYPE.STRING);
//     clusterCol!.name = 'cluster';
//     df.columns.remove('cluster');
//     df.columns.add(clusterCol!);
//     await delay(500);
//   });

//   const viewers = DG.Func.find({package: 'Peptides', tags: ['viewer']}).map((f) => f.friendlyName);
//   for (const v of viewers) {
//     test(v, async () => {
//       await testViewer(v, df.clone(), {detectSemanticTypes: true, arbitraryDfTest: false, readOnly: true});
//       await delay(1000);
//     });
//   }
// });

category('Viewers: Monomer-Position', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  let mpViewer: MonomerPosition;

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    mpViewer = model.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition;

    await delay(500);
  });

  after(async () => await delay(3000));

  test('Tooltip', async () => {
    const cellCoordinates = {col: '9', row: 6};
    const gc = mpViewer.viewerGrid.cell(cellCoordinates.col, cellCoordinates.row);
    const mp = mpViewer.getMonomerPosition(gc);
    expect(showTooltip(model.df, activityCol, Object.entries(model!.settings!.columns!), {
      monomerPosition: mp,
      x: 0,
      y: 0,
      mpStats: model!.monomerPositionStats!,
    }),
    true, `Tooltip is not shown for grid cell at column '${cellCoordinates.col}', row ${cellCoordinates.row}`);
  });

  test('Modes', async () => {
    if (mpViewer === null)
      throw new Error('Monomer-Position viewer doesn\'t exist');

    expect(mpViewer.mode, SELECTION_MODE.MUTATION_CLIFFS,
      `Default Monomer-Position mode is not ${SELECTION_MODE.MUTATION_CLIFFS}`);

    mpViewer.mode = SELECTION_MODE.INVARIANT_MAP;
    expect(mpViewer.mode, SELECTION_MODE.INVARIANT_MAP,
      `Monomer-Position mode is not ${SELECTION_MODE.INVARIANT_MAP} after switching`);

    mpViewer.mode = SELECTION_MODE.MUTATION_CLIFFS;
    expect(mpViewer.mode, SELECTION_MODE.MUTATION_CLIFFS,
      `Monomer-Position mode is not ${SELECTION_MODE.MUTATION_CLIFFS} after switching`);
  });

  test('Highlighting after row removal', async () => {
    const monomerPosition = mpViewer.getMonomerPosition(mpViewer.viewerGrid.cell('9', 6));
    df.rows.removeAt(0, 2);
    await delay(100);

    // Highlighting or's a fresh row-count-sized mask with the cached ones, so stale masks throw 'Array lengths differ'
    highlightMonomerPosition(monomerPosition, df, mpViewer.monomerPositionStats);
  });

  test('Highlighting after row removal with viewer-specific scaling', async () => {
    // Scaling differing from the analysis settings makes the viewer compute stats from its own cached activity column
    mpViewer.getProperty(SAR_PROPERTIES.ACTIVITY_SCALING)!.set(mpViewer, SCALING_METHODS.LG);
    await delay(50);
    expect(mpViewer.monomerPositionStats != null, true, 'Stats should be computed with viewer-specific scaling');
    df.rows.removeAt(0, 2);
    await delay(100);

    const stats = mpViewer.monomerPositionStats;
    const position = Object.keys(stats).find((key) => key !== 'general')!;
    const monomer = Object.keys(stats[position]).find((key) => key !== 'general')!;
    expect(stats[position][monomer].mask.length, df.rowCount,
      'Recomputed masks should match the row count after removal');
    highlightMonomerPosition({monomerOrCluster: monomer, positionOrClusterType: position}, df, stats);
  });
}, {clear: false});

category('Viewers: Most Potent Residues', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  let mprViewer: MostPotentResidues;

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    mprViewer = model.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues;

    await delay(500);
  });

  after(async () => await delay(3000));

  test('Tooltip', async () => {
    const cellCoordinates = {col: 'Diff', row: 6};
    const gc = mprViewer.viewerGrid.cell(cellCoordinates.col, cellCoordinates.row);
    const mp = mprViewer.getMonomerPosition(gc);
    expect(showTooltip(model.df, activityCol, Object.entries(model!.settings!.columns!), {
      monomerPosition: mp,
      x: 0,
      y: 0,
      mpStats: model!.monomerPositionStats!,
    }),
    true, `Tooltip is not shown for grid cell at column '${cellCoordinates.col}', row ${cellCoordinates.row}`);
  });
});

category('Viewers: Logo Summary Table', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  let lstViewer: LogoSummaryTable;

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;

    lstViewer = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable;

    await delay(500);
  });

  after(async () => await delay(3000));

  test('Properties', async () => {
    // Change Logo Summary Table Web Logo Mode property to full
    const webLogoMode = lstViewer.getProperty(LST_PROPERTIES.WEB_LOGO_MODE);
    webLogoMode!.set(lstViewer, PositionHeight.full);
    expect(lstViewer.webLogoMode, PositionHeight.full,
      `Web Logo Mode property is not changed to ${PositionHeight.full}, got ${lstViewer.webLogoMode} instead`);

    // Change Logo Summary Table Members Ratio Threshold proprty to 0
    const threshold = 0;
    const membersRatioThreshold = lstViewer.getProperty(LST_PROPERTIES.MEMBERS_RATIO_THRESHOLD);
    membersRatioThreshold!.set(lstViewer, threshold);
    expect(lstViewer.membersRatioThreshold, threshold,
      `Members Ratio Threshold property is not changed to 0, got ${lstViewer.membersRatioThreshold} instead`);
    expect(lstViewer.viewerGrid.table.filter.anyTrue, true, `Expected to filter out all rows, but ` +
      `${lstViewer.viewerGrid.table.filter.trueCount} rows are left unfiltered`);
  });

  test('Tooltip', async () => {
    const cluster = '0';
    const tooltipElement = lstViewer.showTooltip({
      monomerOrCluster: cluster,
      positionOrClusterType: CLUSTER_TYPE.ORIGINAL,
    }, 0, 0);
    expect(tooltipElement !== null, true, `Tooltip is not shown for cluster '${cluster}'`);
  });

  test('Property change after close', async () => {
    lstViewer.close();
    await delay(100);
    expect(lstViewer.dataFrame == null, true, 'Expected the closed viewer to have no dataFrame');
    lstViewer.getProperty(LST_PROPERTIES.COLUMNS)!.set(lstViewer, [TEST_COLUMN_NAMES.ACTIVITY]);
    lstViewer.render();
  });
}, {clear: false});

category('Viewers: Peptide Generation', () => {
  let df: DG.DataFrame;
  let genViewer: PeptideGenerationViewer;

  before(async () => {
    await PeptideUtils.loadComponents();
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    const sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    genViewer = await df.plot.fromType(PEPTIDE_GENERATION_VIEWER_NAME, {
      [PEPTIDE_GENERATION_PROPS.SEQUENCE]: TEST_COLUMN_NAMES.SEQUENCE,
      [PEPTIDE_GENERATION_PROPS.ACTIVITY]: TEST_COLUMN_NAMES.ACTIVITY,
    }) as PeptideGenerationViewer;
    await delay(500);
  });

  after(async () => await delay(1000));

  test('Generates peptides', async () => {
    const peptideCount = 10;
    genViewer.getProperty(PEPTIDE_GENERATION_PROPS.PEPTIDE_COUNT)!.set(genViewer, peptideCount);
    genViewer.getProperty(PEPTIDE_GENERATION_PROPS.EXCLUDE_EXISTING)!.set(genViewer, false);
    const resultDf = genViewer.generatePeptides();
    expect(resultDf.rowCount > 0, true, 'Peptide Generation produced no peptides');
    expect(resultDf.rowCount <= peptideCount, true,
      `Peptide Generation produced more than the requested ${peptideCount} peptides`);
    for (const colName of Object.values(GEN_COLUMN_NAMES)) {
      expect(resultDf.col(colName) !== null, true,
        `Generated dataframe is missing the '${colName}' column`);
    }
    expect(resultDf.getCol(GEN_COLUMN_NAMES.PEPTIDE).semType, DG.SEMTYPE.MACROMOLECULE,
      'Generated peptide column is not a Macromolecule');
    // Per-position basis columns are appended after the fixed statistics columns.
    expect(resultDf.columns.length > Object.values(GEN_COLUMN_NAMES).length, true,
      'Generated dataframe is missing the per-position basis columns');
  });

  test('Respects High/Low target', async () => {
    genViewer.getProperty(PEPTIDE_GENERATION_PROPS.EXCLUDE_EXISTING)!.set(genViewer, false);
    genViewer.getProperty(PEPTIDE_GENERATION_PROPS.ACTIVITY_TARGET)!.set(genViewer, ACTIVITY_TARGET.HIGH);
    const highDf = genViewer.generatePeptides();
    const highMax = highDf.getCol(GEN_COLUMN_NAMES.PREDICTED).stats.max;

    genViewer.getProperty(PEPTIDE_GENERATION_PROPS.ACTIVITY_TARGET)!.set(genViewer, ACTIVITY_TARGET.LOW);
    const lowDf = genViewer.generatePeptides();
    const lowMin = lowDf.getCol(GEN_COLUMN_NAMES.PREDICTED).stats.min;

    expect(highMax >= lowMin, true,
      'High-target peptides should reach at least as high a predicted activity as the Low-target minimum');
  });
});

category('Viewers: Export', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let split: string[][];

  before(async () => {
    await PeptideUtils.loadComponents();
    df = await loadDemoPeptides();
    split = splitSequences(df.getCol(DEMO_COLUMN_NAMES.SEQUENCE));
    model = await analyzeDemoPeptides(df);
    const mpViewer = model.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition;
    await awaitCheck(() => mpViewer.mutationCliffs !== null, 'Mutation cliffs are not calculated', 30000);
  });

  after(async () => grok.shell.closeAll());

  test('Invariant map', async () => {
    const monomers = [...new Set(split.flat().filter((monomer) => monomer !== ''))].sort();
    const positions = split[0].map((_, i) => `${i + 1}`);
    expect(monomers.length, 22, 'monomers');
    expect(positions.length, 17, 'positions');
    for (const viewerType of [VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP, VIEWER_TYPE.MOST_POTENT_RESIDUES]) {
      (model.findViewer(viewerType) as SARViewer).exportInvariantMap();
      const view = grok.shell.tv;
      const exported = view.dataFrame;
      try {
        expect(exported.name, 'Invariant Map', `${viewerType} export`);
        expect(exported.columns.names().join(', '), [COLUMNS_NAMES.MONOMER, ...positions].join(', '),
          `${viewerType} export columns`);
        expect(exported.getCol(COLUMNS_NAMES.MONOMER).toList().join(', '), monomers.join(', '),
          `${viewerType} export monomers`);
        for (const position of positions) {
          const counts = exported.getCol(position);
          expect(counts.type, DG.COLUMN_TYPE.INT, `${viewerType} export column ${position}`);
          for (let row = 0; row < monomers.length; row++) {
            expect(counts.get(row), split.filter((s) => s[Number(position) - 1] === monomers[row]).length,
              `${viewerType} export count of ${monomers[row]} at ${position}`);
          }
        }
      } finally {
        view.close();
      }
    }
  });

  test('Mutation cliffs', async () => {
    const sequences = df.getCol(DEMO_COLUMN_NAMES.SEQUENCE);
    const activities = df.getCol(DEMO_COLUMN_NAMES.ACTIVITY);
    const ids = df.getCol(DEMO_COLUMN_NAMES.ID);
    const rowByActivity = new Map(activities.toList().map((activity, row) => [activity, row]));
    expect(rowByActivity.size, df.rowCount, 'IC50 identifies the source rows');
    const expectedPairs = singleMutationPairs(split);
    expect(expectedPairs.size, 6253, 'single-mutation pairs');

    (model.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition)['_doExportMutationCliffs']([ids]);
    const view = grok.shell.tv;
    const exported = view.dataFrame;
    try {
      expect(exported.name, 'Mutation Cliffs');
      expect(exported.columns.names().join(', '),
        'Seq 1, Seq 2, Mutation, Seq 1 IC50, Seq 2 IC50, Delta, Seq 1 ID, Seq 2 ID', 'export columns');
      for (const name of ['Seq 1', 'Seq 2']) {
        expect(exported.getCol(name).semType, DG.SEMTYPE.MACROMOLECULE, `${name} semantic type`);
        expect(exported.getCol(name).meta.units, NOTATION.SEPARATOR, `${name} units`);
      }
      expect(exported.getCol('Mutation').semType, SEM_TYPES.MACROMOLECULE_DIFFERENCE, 'Mutation semantic type');
      const exportedPairs = new Set<string>();
      for (let row = 0; row < exported.rowCount; row++) {
        const activity1 = exported.get('Seq 1 IC50', row);
        const activity2 = exported.get('Seq 2 IC50', row);
        const a = rowByActivity.get(activity1)!;
        const b = rowByActivity.get(activity2)!;
        const pair = `${Math.min(a, b)}:${Math.max(a, b)}`;
        expect(expectedPairs.has(pair) && !exportedPairs.has(pair), true, `row ${row}: duplicate or non-cliff pair`);
        exportedPairs.add(pair);
        const seq1 = exported.get('Seq 1', row);
        const seq2 = exported.get('Seq 2', row);
        expect(seq1, sequences.get(a), `row ${row}: Seq 1`);
        expect(seq2, sequences.get(b), `row ${row}: Seq 2`);
        expect(exported.get('Mutation', row), `${seq1}#${seq2}`, `row ${row}: Mutation`);
        expectFloat(exported.get('Delta', row), activity1 - activity2, 1e-12, `row ${row}: Delta`);
        expect(exported.get('Seq 1 ID', row), ids.get(a), `row ${row}: Seq 1 ID`);
        expect(exported.get('Seq 2 ID', row), ids.get(b), `row ${row}: Seq 2 ID`);
      }
      expect(exportedPairs.size, expectedPairs.size, 'exported single-mutation pairs');
    } finally {
      view.close();
    }
  });
}, {clear: false});

category('Viewers: Mutation Cliffs', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let split: string[][];

  before(async () => {
    await PeptideUtils.loadComponents();
    df = await loadDemoPeptides(200);
    split = splitSequences(df.getCol(DEMO_COLUMN_NAMES.SEQUENCE));
    model = await analyzeDemoPeptides(df);
  });

  after(async () => grok.shell.closeAll());

  test('Cliff counts', async () => {
    const mpViewer = model.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition;
    await awaitCheck(() => mpViewer.mutationCliffs !== null, 'Mutation cliffs are not calculated', 30000);
    const cliffs = mpViewer.mutationCliffs!;
    const cliffCount = (monomer: string, position: string): number =>
      Array.from(cliffs.get(monomer)?.get(position)?.values() ?? []).reduce((sum, targets) => sum + targets.length, 0);
    let pairs = 0;
    const uniquePairs = new Set<string>();
    for (const [monomer, positions] of cliffs) {
      for (const [position, indexes] of positions) {
        pairs += cliffCount(monomer, position);
        for (const [from, targets] of indexes) {
          for (const to of targets)
            uniquePairs.add(`${Math.min(from, to)}:${Math.max(from, to)}`);
        }
      }
    }
    expect(pairs, 4484, 'cliff pairs, one per mutated monomer of a pair');
    expect(uniquePairs.size, 2242, 'unique cliff pairs');
    const expectedPairs = singleMutationPairs(split);
    expect(expectedPairs.size === uniquePairs.size && [...expectedPairs].every((pair) => uniquePairs.has(pair)), true,
      'cliff pairs are the single-mutation pairs');
    expect(cliffCount('A', '2'), 100, 'cliffs of A at 2');
    expect(cliffCount('A', '3'), 0, 'cliffs of A at 3');

    expect(model.positionColumns?.length, 17, 'positions');
    expect(df.getCol('2').semType, SEM_TYPES.MONOMER, 'position column semantic type');
    const stats = mpViewer.monomerPositionStats['2']['A'];
    expect(stats.count, 59, 'count of A at 2');
    expect(Number.isFinite(stats.meanDifference), true, 'mean difference of A at 2');
    expect(stats.pValue !== null && stats.pValue >= 0 && stats.pValue <= 1, true, `p-value ${stats.pValue} of A at 2`);
  });

  test('Position chart', async () => {
    const participants = (position: number): number[] => {
      const rows = new Set<number>();
      for (const pair of singleMutationPairs(split)) {
        const [a, b] = pair.split(':').map(Number);
        if (split[a][position - 1] !== split[b][position - 1])
          rows.add(a).add(b);
      }
      return [...rows].sort((a, b) => a - b);
    };
    const viewer = await df.plot.fromType(VIEWER_TYPE.SEQUENCE_MUTATION_CLIFFS, {
      sequenceColumnName: DEMO_COLUMN_NAMES.SEQUENCE, activityColumnName: DEMO_COLUMN_NAMES.ACTIVITY, position: 2,
    }) as MutationCliffsViewer;
    model.analysisView.addViewer(viewer);
    const chartAt = async (position: number): Promise<DG.DataFrame> => {
      viewer.getProperty('position')!.set(viewer, position);
      await awaitCheck(() => viewer.getWidgetStatus().values.position === position && !viewer.isRenderPending,
        `The position ${position} chart is not rendered`, 15000);
      return viewer.innerDf;
    };

    const rows = participants(2);
    expect(rows.length, 115, 'peptides in position 2 cliffs');
    const chart = await chartAt(2);
    expect(chart.rowCount, rows.length, 'position 2 chart rows');
    const monomers = chart.getCol('Position 2');
    const plotted = chart.getCol(viewer.activityColumnName);
    const activities = df.getCol(DEMO_COLUMN_NAMES.ACTIVITY);
    for (let i = 0; i < rows.length; i++) {
      expect(monomers.get(i) ?? '', split[rows[i]][1], `position 2 chart row ${i}: monomer`);
      expect(Math.fround(plotted.get(i)), Math.fround(activities.get(rows[i])), `position 2 chart row ${i}: activity`);
    }
    expect(participants(1).length, 0, 'peptides in position 1 cliffs');
    expect((await chartAt(1)).rowCount, 0, 'position 1 chart rows');
  });
}, {clear: false});
