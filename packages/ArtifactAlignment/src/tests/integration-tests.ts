import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

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
