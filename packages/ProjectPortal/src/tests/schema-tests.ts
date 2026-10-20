import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import dayjs from 'dayjs';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, stagesOf} from '../portal';
import {_package} from '../package-test';

/** Schema 0.2: the run-selection loop and the workspace vocabulary, round-tripped through the
 * domain API on real rows. Every fixture carries the TEST- prefix and is removed afterwards. */
category('Portal: runs and criteria', () => {
  test('a saved run is selected for a document criterion', async () => {
    await loadAllPortfolios(_package);
    const cmc = await portfolioByCode('EX-BIOC');
    const stage = (await stagesOf(cmc.stage_model_id))[0];
    const tag = `TEST-${Date.now() % 1e9}`;
    const undo: (() => Promise<unknown>)[] = [];
    try {
      const abs = DG.Func.find({name: 'Abs'})[0];
      const call = abs.prepare({[abs.inputs[0].name]: -2});
      await call.call();
      const saved = await grok.dapi.functions.calls.save(call);
      undo.push(() => grok.dapi.functions.calls.delete(saved));

      const [program] = await portalDb.programs.insert({portfolio_id: cmc.id, code: tag, name: tag,
        stage_model_id: cmc.stage_model_id, stage_id: stage.id, source: 'manual', modality: 'MAB'});
      undo.push(() => portalDb.programs.delete(program.id));
      const [op] = await portalDb.unitOperations.insert({code: `${tag}-OP`, name: 'Test operation'});
      undo.push(() => portalDb.unitOperations.delete(op.id));
      const [runType] = await portalDb.runTypes.insert({code: `${tag}-RT`, name: 'Simulated'});
      undo.push(() => portalDb.runTypes.delete(runType.id));
      const [model] = await portalDb.models.insert({code: `${tag}:Abs`, name: 'Test model'});
      undo.push(() => portalDb.models.delete(model.id));
      const [criterion] = await portalDb.gateCriterions.insert({stage_id: stage.id, code: `${tag}-DOC`,
        name: 'Test deliverable', kind: 'document', modality: 'MAB'});
      undo.push(() => portalDb.gateCriterions.delete(criterion.id));

      const [status] = await portalDb.criterionStatuses.insert({program_id: program.id,
        gate_criterion_id: criterion.id, status: 'not_started'});
      const [run] = await portalDb.modelRuns.insert({func_call: saved.id, program_id: program.id,
        model_id: model.id, model_name: abs.nqName, title: `${tag} run`, stage_id: stage.id,
        unit_operation_id: op.id, run_type_id: runType.id, purpose: 'characterization'});
      undo.push(() => portalDb.modelRuns.delete(run.id));

      const [pick] = await portalDb.criterionArtifacts.insert({program_id: program.id,
        criterion_status_id: status.id, model_run_id: run.id, role: 'candidate'});
      await portalDb.criterionArtifacts.update(pick.id, {role: 'selected', selected_by: grok.shell.user.id,
        selected_on: dayjs(), rationale: 'closest to the golden batch'});
      await portalDb.criterionStatuses.update(status.id, {status: 'met'});

      const picks = await portalDb.criterionArtifacts.query({filter: `criterion_status_id = "${status.id}"`});
      expect(picks.length, 1);
      expect(picks[0].role, 'selected');
      const met = await portalDb.criterionStatuses.query({filter: `id = "${status.id}"`});
      expect(met[0].status, 'met');
      const runs = await portalDb.modelRuns.query({filter: `program_id = "${program.id}"`});
      expect(runs.length, 1);
      expect(runs[0].func_call, saved.id);

      let refused = false;
      try {
        await portalDb.modelRuns.update(run.id, {stage_id: (await stagesOf(cmc.stage_model_id))[1].id});
      } catch {
        refused = true;
      }
      expect(refused, true, 'the stage at run is write-once');
    } finally {
      for (const u of undo.reverse())
        await u().catch(() => {});
    }
  });

  test('a program runs one flow per process type', async () => {
    await loadAllPortfolios(_package);
    const cmc = await portfolioByCode('EX-BIOC');
    const tag = `TEST-${Date.now() % 1e9}`;
    const undo: (() => Promise<unknown>)[] = [];
    try {
      const [op] = await portalDb.unitOperations.insert({code: `${tag}-THAW`, name: 'Thaw'});
      undo.push(() => portalDb.unitOperations.delete(op.id));
      const [api] = await portalDb.processFlows.insert({code: `${tag}-API`, name: 'Test API', process_type: 'Drug Substance'});
      undo.push(() => portalDb.processFlows.delete(api.id));
      await portalDb.processSteps.insert({process_flow_id: api.id, unit_operation_id: op.id, step_number: 1,
        folder_name: '01 - Thaw'});
      const [rule] = await portalDb.flowRules.insert({code: `MAB:${tag}-API`, modality: 'MAB',
        process_flow_id: api.id, priority: 1});
      undo.push(() => portalDb.flowRules.delete(rule.id));
      const [program] = await portalDb.programs.insert({portfolio_id: cmc.id, code: tag, name: tag,
        stage_model_id: cmc.stage_model_id, source: 'manual', modality: 'MAB'});
      undo.push(() => portalDb.programs.delete(program.id));

      await portalDb.programFlows.insert({program_id: program.id, process_flow_id: api.id, source: 'rule'});
      const flows = await portalDb.programFlows.query({filter: `program_id = "${program.id}"`});
      expect(flows.length, 1);
      expect(flows[0].source, 'rule');
      const steps = await portalDb.processSteps.query({filter: `process_flow_id = "${api.id}"`});
      expect(steps.map((s) => s.folder_name).join(), '01 - Thaw');
    } finally {
      for (const u of undo.reverse())
        await u().catch(() => {});
    }
  });
});
