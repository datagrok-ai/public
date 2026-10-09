import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {getToolTip, TreeNode} from '../components/TreeWizard/TreeNode';

category('TreeWizard: step tooltip', () => {
  test('Failed step tooltip shows the run error', async () => {
    expect(getToolTip('failed', false, undefined, undefined, undefined, 'Error: division by zero'),
      'Run failed: Error: division by zero');
  });

  test('Failed step without a stored error keeps the plain text', async () => {
    expect(getToolTip('failed', false), 'Run failed');
  });
});

category('TreeWizard: step node', () => {
  const mountNode = (isSelected: boolean) => {
    const root = ui.div();
    document.body.append(root);
    const app = Vue.createApp(TreeNode, {
      stat: {data: {uuid: 'u1', type: 'funccall', configId: 'step1'}, children: [], parent: null, open: false} as any,
      callState: {isRunning: false, isRunnable: true, isOutputOutdated: true, pendingDependencies: []} as any,
      isDeletable: true,
      isSelected,
    });
    app.mount(root);
    return {root, unmount: () => {
      app.unmount();
      root.remove();
    }};
  };

  test('The selected node shows its actions with labels', async () => {
    const selected = mountNode(true);
    const other = mountNode(false);
    try {
      expect(selected.root.querySelector('[aria-label="Remove"]')?.getAttribute('role'), 'button');
      expect(other.root.querySelector('[aria-label="Remove"]'), null);
      expect(selected.root.querySelector('[aria-label="This step is available to run"]')?.getAttribute('role'), 'img');
    } finally {
      selected.unmount();
      other.unmount();
    }
  });
});
