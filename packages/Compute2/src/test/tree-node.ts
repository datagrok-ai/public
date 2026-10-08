import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {getToolTip} from '../components/TreeWizard/TreeNode';

category('TreeWizard: step tooltip', () => {
  test('Failed step tooltip shows the run error', async () => {
    expect(getToolTip('failed', false, undefined, undefined, undefined, 'Error: division by zero'),
      'Run failed: Error: division by zero');
  });

  test('Failed step without a stored error keeps the plain text', async () => {
    expect(getToolTip('failed', false), 'Run failed');
  });
});
