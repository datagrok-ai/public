import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {T_ALIGNMENT, T_HISTORY} from '../domain/constants';
import {driftCheck} from '../service/drift-check';
import {getLiveVersions, publishWorkflowRun} from '../service/publication-service';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram} from './fixtures';

category('ArtifactAlignment: drift check', () => {
  test('healthy catalog reports no findings for its publications', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Healthy'});
    await publishWorkflowRun({sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Healthy'});
    const findings = await driftCheck();
    expect(findings.filter((f) => f.publicationId === v1.publicationId).length, 0,
      JSON.stringify(findings));
  });

  test('detects a publication left with no live row (interrupted republish)', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'NoLive'});
    const row: any = (await getLiveVersions(v1.publicationId))[0];
    await grok.dapi.domains.table(T_HISTORY).insert({
      publication_id: row.publication_id, revision: row.revision, name: row.name,
      artifact_id: row.artifact_id, program_id: row.program_id, final_status: row.status,
    });
    await grok.dapi.domains.table(T_ALIGNMENT).delete(row.id);
    const findings = await driftCheck();
    expect(findings.some((f) => f.kind === 'no-live-row' && f.publicationId === v1.publicationId),
      true, JSON.stringify(findings));
  });

  test('detects an interrupted approve flip (in-review live row, approved only in history)', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Flip'});
    const v2 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Flip', skipAutoApprove: true});
    // Simulate the crash between archive and flip: archive v1 without approving v2.
    const approved: any = (await getLiveVersions(v1.publicationId)).find((r) => r.status === 'approved');
    await grok.dapi.domains.table(T_HISTORY).insert({
      publication_id: approved.publication_id, revision: approved.revision, name: approved.name,
      artifact_id: approved.artifact_id, program_id: approved.program_id, final_status: 'approved',
    });
    await grok.dapi.domains.table(T_ALIGNMENT).delete(approved.id);
    const findings = await driftCheck();
    expect(findings.some((f) => f.kind === 'interrupted-flip' && f.publicationId === v2.publicationId),
      true, JSON.stringify(findings));
  });

  test('deep mode verifies the frozen artifact funccalls resolve', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Deep'});
    const findings = await driftCheck({deep: true});
    expect(findings.some((f) => f.kind === 'missing-frozen-artifact' &&
      f.publicationId === v1.publicationId), false, JSON.stringify(findings));
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});

category('ArtifactAlignment: integration surface', () => {
  test('publish functions are registered for the Compute2 gate to find', async () => {
    expect(DG.Func.find({name: 'publishWorkflowRunDialog'}).length > 0, true);
    expect(DG.Func.find({name: 'publishWorkflowRun'}).length > 0, true);
  });

  test('Compute2 exposes the open-by-id function used by the catalog', async () => {
    // Present only when the matching Compute2 build is deployed alongside.
    expect(DG.Func.find({name: 'OpenWorkflowRun'}).length > 0, true,
      'Compute2 with OpenWorkflowRun is not deployed');
  });
});
