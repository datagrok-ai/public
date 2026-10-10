import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after as afterAll, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programContext, programData} from '../portal';
import {
  compareVersions, compoundColumn, describeFit, describePrediction, latestPredictions, planScoring, Prediction,
  predictedVsMeasured, predictionColumn, unitsOf, unscored, withPredictions,
} from '../model/predictions';
import {predictionsOf, scoreCompounds} from '../predictions';
import {recommendFor} from '../recommend';
import {provenanceOf} from '../provenance';
import {deleteRunTable} from '../model-runs';
import {_package} from '../package-test';

const SOL = 'ProjectPortalCore:SolubilityEstimate';

const p = (compound: string, value: number | null, version = '1.0', on = '2026-09-01T00:00:00Z'): Prediction =>
  ({compound, predicts: 'sol_um', model: 'M', modelName: 'Model M', version, value, on});

category('Portal: prediction rules', () => {
  test('the output columns: the estimate for a column, never the measured copy; the compound; units by suffix', async () => {
    const cols = ['compound', 'solubility_um', 'estimated_solubility_um'];
    expect(predictionColumn(cols, 'solubility_um'), 'estimated_solubility_um');
    expect(predictionColumn(['Compound', 'herg_ic50_um_pred'], 'herg_ic50_um'), 'herg_ic50_um_pred');
    expect(predictionColumn(['compound', 'solubility_um'], 'solubility_um'), null);
    expect(compoundColumn(['row', 'Compound_ID']), 'Compound_ID');
    expect(compoundColumn(['x']), null);
    expect(`${unitsOf('solubility_um')} ${unitsOf('ic50_nM')} ${unitsOf('papp')}`, 'uM nM null');
  });

  test('scoring stores one value per compound the version has not scored; again, nothing', async () => {
    const scored = [{compound: 'C1', value: 10}, {compound: 'C2', value: '12.5'}, {compound: null, value: 3},
      {compound: 'C3', value: null}, {compound: 'C1', value: 11}];
    const first = planScoring(scored, [], 'M', '1.0', 'sol_um');
    expect(first.insert.map((x) => `${x.compound}=${x.value}`).join(), 'C1=10,C2=12.5');
    expect(`${first.kept} ${first.skipped}`, '0 2');
    const again = planScoring(scored, [p('C1', 10), p('C2', 12.5)], 'M', '1.0', 'sol_um');
    expect(`${again.insert.length} ${again.kept}`, '0 2');
    expect(planScoring(scored, [p('C1', 10)], 'M', '2.0', 'sol_um').insert.length, 2, 'a new version scores again');
    expect(unscored(['C1', 'C2', 'C3'], [p('C1', 1)], 'M', '1.0', 'sol_um').join(), 'C2,C3');
  });

  test('the latest prediction per compound: the newest version, numbers compared as numbers', async () => {
    expect(compareVersions('1.10', '1.9') > 0 && compareVersions('2.0', '2.0') === 0 && compareVersions('1.0', '1.0.1') < 0, true);
    const latest = latestPredictions([p('C1', 1, '1.9'), p('C1', 2, '1.10'), p('C2', 3, '1.0'), p('C2', 4, '1.0', '2026-09-02T00:00:00Z')]);
    expect(latest.map((x) => `${x.compound}=${x.value}`).join(), 'C1=2,C2=4');
    expect(describePrediction({...p('C1', 34.0999), predicts: 'solubility_um'}), '34.1 uM (Model M 1.0)');
  });

  test('predicted against measured: pairs, mean absolute error, r; missing assays take the stored predictions', async () => {
    const rows = [{compound: 'C1', sol_um: 10}, {compound: 'C2', sol_um: 20}, {compound: 'C3', sol_um: 30}, {compound: 'C4', sol_um: null}];
    const [f] = predictedVsMeasured([p('C1', 12), p('C2', 18), p('C3', 33), p('C4', 25)], rows);
    expect(`${f.pairs.length} ${f.paired} ${f.mae!.toFixed(2)} ${f.r!.toFixed(3)}`, '4 3 2.33 0.971');
    expect(describeFit(f), '3 of 4 measured; mean absolute error 2.3 uM; r 0.97');
    const missing = [{criterion: 'adme', criterionName: 'ADME panel', required: true, column: 'sol_um', missing: 1, items: ['C4']}];
    const [m] = withPredictions(missing, [p('C1', 12), p('C4', 25)]);
    expect(m.predicted.map((x) => `${x.compound}=${x.value}`).join(), 'C4=25');
  });
});

category('Portal: predictions', () => {
  let smd: PortfolioRow;
  let kin1: ProgramRow;
  let first: Awaited<ReturnType<typeof scoreCompounds>> | null = null;

  /** Removes the example model's predictions on EXS-KIN1, and the scoring runs behind them. */
  const clear = async () => {
    const rows = await portalDb.predictions.query({filter: `program_id = "${kin1.id}" and model_code = "${SOL}"`});
    const runs = [...new Set(rows.map((r) => r.model_run_id).filter((x): x is string => !!x))];
    for (const r of rows)
      await portalDb.predictions.delete(r.id).catch(() => undefined);
    for (const id of runs) {
      const run = await portalDb.modelRuns.get(id).catch(() => null);
      if (!run)
        continue;
      const call = run.run_id ? await grok.dapi.functions.calls.find(run.run_id).catch(() => null) : null;
      for (const t of (call?.options?.['portalTables'] as string[] | undefined) ?? [])
        await deleteRunTable(t).catch(() => undefined);
      if (call)
        await grok.dapi.functions.calls.delete(call).catch(() => undefined);
      await portalDb.modelRuns.delete(id).catch(() => undefined);
    }
  };

  before(async () => {
    await loadAllPortfolios(_package);
    smd = await portfolioByCode('EX-SMD');
    kin1 = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-KIN1"`}))[0];
    await clear();
  });

  afterAll(async () => {
    await clear();
  });

  test('a scoring runs once and stores one prediction per compound, each pointing at the saved run', async () => {
    first = await scoreCompounds('EX-SMD', 'EXS-KIN1', SOL);
    expect(first.stored.map((s) => `${s.predicts}=${s.count}`).join(), 'solubility_um=11', first.summary);
    expect(first.version, '1.0');
    const {rows} = await predictionsOf(kin1.id);
    const mine = rows.filter((r) => r.model_code === SOL);
    expect(mine.length, 11);
    expect(mine.every((r) => r.model_run_id === first!.modelRunId && r.units === 'uM' && r.model_version === '1.0'), true);
    const run = await portalDb.modelRuns.get(first.modelRunId!);
    expect(run.model_name, SOL);
    expect(run.run_id, first.runId);
  });

  test('scoring the same version again runs nothing and stores nothing', async () => {
    const again: string = await grok.functions.call('ProjectPortalCore:scoreCompounds',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', model: SOL});
    expect(again, 'EXS-KIN1: all 11 compounds already scored by Solubility estimate (example) 1.0; nothing run');
    expect((await portalDb.predictions.query({filter: `program_id = "${kin1.id}" and model_code = "${SOL}"`})).length, 11);
  });

  test('predicted vs measured on solubility_um; the missing assay shows the stored prediction', async () => {
    const {latest} = await predictionsOf(kin1.id);
    const rows = (await programData(smd, 'EXS-KIN1'))!;
    const fit = predictedVsMeasured(latest.filter((x) => x.model === SOL), rows)[0];
    expect(`${fit.pairs.length} ${fit.paired}`, '11 7');
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:programPredictions',
      {portfolio: 'EX-SMD', program: 'EXS-KIN1', model: SOL});
    expect(df.rowCount, 11);
    expect(df.col('measured')!.stats.valueCount, 7);
    const ctx = await programContext(kin1.id);
    const lo = ctx.stages.find((s) => s.code === 'LO')!;
    const r = await recommendFor(ctx, {stageId: lo.id});
    const sol = r.missing.find((m) => m.column === 'solubility_um')!;
    const [m] = withPredictions([sol], latest);
    expect(m.predicted.map((x) => `${x.compound} ${describePrediction(x)}`).join(), 'EXC-1006 34.1 uM (Solubility estimate (example) 1.0)');
  });

  test('provenance leads from the scoring run to the model and the predictions', async () => {
    // the first test's run, else the run behind the stored predictions: this test stands on its own
    const runId = first?.modelRunId ?? (await portalDb.predictions.query({filter: `program_id = "${kin1.id}" and model_code = "${SOL}"`,
      limit: 1}))[0]?.model_run_id;
    expect(!!runId, true, 'no scoring run of the example model on EXS-KIN1');
    const prov = await provenanceOf('model_run', runId!, 2);
    const up = prov.nodes.filter((n) => n.direction === 'up').map((n) => n.kind);
    const down = prov.nodes.filter((n) => n.direction === 'down').map((n) => `${n.kind}:${n.name}`);
    expect(up.includes('model'), true, up.join());
    expect(down.includes('prediction:solubility_um on 11 compounds'), true, down.join());
  });
});
