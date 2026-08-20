import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {driftCheck} from '../service/drift-check';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram, publishRun} from './fixtures';

category('ArtifactAlignment: drift check', () => {
  test('healthy catalog reports no findings for its publications', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Healthy'});
    await publishRun({sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Healthy'});
    const findings = await driftCheck();
    expect(findings.filter((f) => f.publicationId === v1.publicationId).length, 0,
      JSON.stringify(findings));
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
