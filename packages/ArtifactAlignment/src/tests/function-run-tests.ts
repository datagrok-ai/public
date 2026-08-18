import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {cloneRun} from '../service/clone-service';
import {
  OPT_FROZEN, OPT_PARENT_CALL_ID, OPT_PUBLICATION_ID, T_ALIGNMENT,
} from '../domain/constants';
import {discover, getPublicationHistory, publishWorkflowRun} from '../service/publication-service';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram} from './fixtures';

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

  test('the frozen function run is loadable, stamped, and audience-granted', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun(7);
    const result = await publishWorkflowRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'Frozen single'});
    const frozen = await historyUtils.loadRun(result.artifactId);
    expect(frozen.options[OPT_FROZEN], 'true');
    expect(frozen.options[OPT_PUBLICATION_ID], result.publicationId);
    expect(frozen.options[OPT_PARENT_CALL_ID] == null, true, 'a single run has no parent call');
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

  test('frozen copies stay out of the author history listings', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const workflow = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Hidden wf'});
    const single = await publishWorkflowRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'Hidden single'});
    const author = await grok.dapi.users.current();
    const providerRuns = await historyUtils.pullRunsByName(
      'AaTestProvider', [{author: author as any}], {}, ['options']);
    expect(providerRuns.some((r) => r.id === workflow.artifactId), false,
      'the frozen meta call must not appear in the workflow history list');
    expect(providerRuns.some((r) => r.id === run.metaCallId), true,
      'the source run must still appear');
    const stepRuns = await historyUtils.pullRunsByName(
      'AaTestStep', [{author: author as any}], {}, ['options']);
    expect(stepRuns.some((r) => r.id === single.artifactId), false,
      'the frozen single run must not appear in the step history list');
  });

  test('cloneRun detects the run kind', async () => {
    const run = await makeSavedRun();
    const single = await cloneRun(run.stepCallId, {publicationId: crypto.randomUUID()});
    expect(single.artifactType, 'function-run');
    expect(Object.keys(single.childIdMap).length, 0);
    const workflow = await cloneRun(run.metaCallId, {publicationId: crypto.randomUUID()});
    expect(workflow.artifactType, 'workflow-run');
    expect(Object.keys(workflow.childIdMap).length, 1);
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
