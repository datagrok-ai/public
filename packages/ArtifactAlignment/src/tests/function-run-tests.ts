import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {cloneRun} from '../service/clone-service';
import {
  OPT_FROZEN, OPT_PUBLICATION_ID, T_ALIGNMENT,
} from '../domain/constants';
import {discover, getPublicationHistory, publishWorkflowRun} from '../service/publication-service';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram, STEP_NQ} from './fixtures';

// Single saved runs — RFV model runs and individually saved workflow steps — publish
// through the same flow as workflow runs, detected by the absence of a pipeline config.
category('ArtifactAlignment: function runs', () => {
  test('a single saved run publishes as a function-run artifact', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun(4);
    const result = await publishWorkflowRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'Step artifact',
      workstream: 'modeling'});
    expect(result.status, 'approved');
    expect(result.artifactId === run.stepCallId, false, 'the frozen copy must be a new call');
    const discovered = await discover(program.id);
    expect(discovered.length, 1);
    expect(discovered[0].artifact_type, 'function-run');
  });

  test('a live in-memory run publishes without a prior history save', async () => {
    const program = await makeTestProgram();
    const fc = DG.Func.byName(STEP_NQ).prepare({a: 5});
    await fc.call(false, undefined, {processed: true, report: false});
    const result = await publishWorkflowRun({
      sourceCall: fc, programId: program.id, name: 'In-memory artifact'});
    expect(result.status, 'approved');
    const frozen = await historyUtils.loadRun(result.artifactId);
    expect(frozen.options[OPT_FROZEN], 'true');
    expect((frozen.outputs['result'] as DG.DataFrame).rowCount, 5);
    expect(fc.options[OPT_FROZEN] == null, true, 'the live call must not be mutated');
  });

  test('the frozen function run is loadable, stamped, and audience-granted', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun(7);
    const result = await publishWorkflowRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'Frozen single'});
    const frozen = await historyUtils.loadRun(result.artifactId);
    expect(frozen.options[OPT_FROZEN], 'true');
    expect(frozen.options[OPT_PUBLICATION_ID], result.publicationId);
    expect((frozen.outputs['result'] as DG.DataFrame).rowCount, 7);
    const perms: any = await grok.dapi.permissions.get(
      (frozen.outputs['result'] as DG.DataFrame).getTableInfo());
    const grantedIds: string[] = [];
    for (const key of Object.keys(perms ?? {})) {
      for (const g of perms[key] ?? [])
        grantedIds.push(g?.id);
    }
    expect(grantedIds.includes(program.viewers.id), true, 'viewers group must hold View on the dataframe');
  });

  test('republishing a function run under the same key bumps the revision', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'Single repub'});
    const v2 = await publishWorkflowRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'Single repub'});
    expect(v1.publicationId, v2.publicationId);
    expect(v2.revision, 2);
    expect((await getPublicationHistory(v1.publicationId))[0].revision, 1);
  });

  test('workflow and function runs coexist in discovery with distinct type facets', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'A workflow'});
    await publishWorkflowRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'A step'});
    const facets: any = await grok.dapi.domains.table(T_ALIGNMENT).facets({
      filter: [{property: 'program_id', operator: '=', value: program.id}],
      facets: [{id: 'type', kind: 'categories', column: 'artifact_type'}],
    });
    const counts = Object.fromEntries(
      facets.facets.type.categories.map((c: any) => [c.value, c.filtered]));
    expect(counts['workflow-run'], 1, JSON.stringify(counts));
    expect(counts['function-run'], 1, JSON.stringify(counts));
  });

  test('frozen copies are excluded by the Compute2 history frozen filter', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const workflow = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Hidden wf'});
    const single = await publishWorkflowRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'Hidden single'});
    // The persisted options drive the listing filter in Compute2's History component;
    // fetch the involved runs by id — an unbounded history pull is not what's under test.
    const fetch = (id: string) =>
      grok.dapi.functions.calls.allPackageVersions().include('options').find(id);
    const runs = await Promise.all(
      [workflow.artifactId, single.artifactId, run.metaCallId].map(fetch));
    const visible = runs.filter((r) => r.options[OPT_FROZEN] !== 'true');
    expect(visible.length, 1, 'both frozen copies must be filtered out');
    expect(visible[0].id, run.metaCallId, 'the source run must pass the filter');
  });

  test('cloneRun detects the run kind', async () => {
    const run = await makeSavedRun();
    const single = await cloneRun(run.stepCallId, {publicationId: crypto.randomUUID()});
    expect(single.artifactType, 'function-run');
    const workflow = await cloneRun(run.metaCallId, {publicationId: crypto.randomUUID()});
    expect(workflow.artifactType, 'workflow-run');
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
