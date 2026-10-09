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

  test('A locked step says so instead of next', async () => {
    expect(getToolTip('next', true), 'This step is locked');
    expect(getToolTip('succeeded', true), 'This step is succeeded');
  });

  test('Tooltips list the inputs behind the status', async () => {
    const validations = {a: {errors: ['bad']}, b: {warnings: ['careful']}};
    const consistency = {
      c: {restriction: 'restricted', inconsistent: true, assignedValue: 1},
      d: {restriction: 'info', inconsistent: true, assignedValue: 1},
    } as const;
    expect(getToolTip('next error', false, undefined, validations, consistency), 'This step needs user input: a');
    expect(getToolTip('next warn', false, undefined, validations, consistency),
      'This step is available to run, but has warnings: b, c');
    expect(getToolTip('succeeded warn', false, undefined, validations, consistency),
      'This step is succeeded, but has warnings: a, b');
    expect(getToolTip('succeeded info', false, undefined, validations, consistency),
      'This step is succeeded with changes: d');
    expect(getToolTip('succeeded inconsistent', false, undefined, validations, consistency),
      'This step is succeeded, but has inconsistent inputs: c, d');
    expect(getToolTip('succeeded', false, undefined, validations, consistency), 'This step is succeeded');
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
