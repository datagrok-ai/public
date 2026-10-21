import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {getToolTip, TreeNode} from '../components/TreeWizard/TreeNode';
import {sourceRunTooltip} from '../components/TreeWizard/source-run';

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

  test('A loaded node shows one source run icon: a lock on the left if read-only, else last in the row', async () => {
    const mountNode = (data: Record<string, any>, isReadonly = false) => {
      const root = ui.div();
      document.body.append(root);
      const app = Vue.createApp(TreeNode, {
        stat: {data: {uuid: 'n1', configId: 'n', ...data}, children: [], parent: null, open: false} as any,
        callState: data.type === 'funccall' ?
          {isRunning: false, isRunnable: true, isOutputOutdated: false, pendingDependencies: []} as any : undefined,
        isReadonly,
        isDeletable: true,
        isSelected: true,
      });
      app.mount(root);
      return {root, unmount: () => {
        app.unmount();
        root.remove();
      }};
    };
    const basedOn = mountNode({type: 'static', steps: [], sourceRunId: 'run1'});
    const basedOnStep = mountNode({type: 'funccall', sourceRunId: 'run4', sourceRunTitle: 'Baseline'});
    const asIs = mountNode({type: 'funccall', isReadonly: true, sourceRunId: 'run2'}, true);
    const asIsPipeline = mountNode({type: 'static', steps: [], isReadonly: true, sourceRunId: 'run3'}, true);
    const plain = mountNode({type: 'static', steps: []});
    const icon = (node: {root: HTMLElement}, isReadonly: boolean) =>
      node.root.querySelector(`dg-icon-fa[aria-label*="${sourceRunTooltip(isReadonly)}"]`);
    try {
      expect(!!icon(basedOn, false)?.querySelector('.fa-code-branch'), true);
      expect(!!basedOn.root.querySelector('dg-icon-fa[aria-label="Remove"]'), true);
      expect(icon(basedOn, false)?.nextElementSibling, null);
      expect(!!basedOnStep.root.querySelector(`dg-icon-fa[aria-label="${sourceRunTooltip(false, 'Baseline')}"]`), true);
      for (const node of [asIs, asIsPipeline]) {
        expect(!!icon(node, true)?.querySelector('.fa-lock'), true);
        expect(!!icon(node, true)?.closest('.ml-auto'), false);
        expect(node.root.querySelectorAll('dg-icon-fa .fa-lock').length, 1);
        expect(node.root.querySelector('.fa-code-branch'), null);
      }
      expect(plain.root.querySelector('.fa-code-branch, .fa-lock'), null);
    } finally {
      basedOn.unmount();
      basedOnStep.unmount();
      asIs.unmount();
      asIsPipeline.unmount();
      plain.unmount();
    }
  });

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
