import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import {category, test, expect, expectArray, waitUntil} from '@datagrok-libraries/test/src/test';
import {getToolTip, TreeNode} from '../components/TreeWizard/TreeNode';
import {sourceRunLabel, sourceRunTooltip} from '../components/TreeWizard/source-run';

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
  const mountNode = (
    data: Record<string, any>, {isReadonly = false, isSelected = true, isOutputOutdated = false} = {},
  ) => {
    const root = ui.div();
    document.body.append(root);
    const app = Vue.createApp(TreeNode, {
      stat: {data: {uuid: 'n1', configId: 'n', ...data}, children: [], parent: null, open: false} as any,
      callState: data.type === 'funccall' ?
        {isRunning: false, isRunnable: true, isOutputOutdated, pendingDependencies: []} as any : undefined,
      isReadonly,
      isDeletable: true,
      isSelected,
    });
    app.mount(root);
    return {root, unmount: () => {
      app.unmount();
      root.remove();
    }};
  };

  test('A loaded node shows the source run menu icon last in the row, a read-only one keeps its lock', async () => {
    const basedOn = mountNode({type: 'static', steps: [], sourceRunId: 'run1'});
    const basedOnStep = mountNode({type: 'funccall', sourceRunId: 'run4', sourceRunTitle: 'Baseline'});
    const asIs = mountNode({type: 'funccall', isReadonly: true, sourceRunId: 'run2'}, {isReadonly: true});
    const asIsPipeline =
      mountNode({type: 'static', steps: [], isReadonly: true, sourceRunId: 'run3'}, {isReadonly: true});
    const plain = mountNode({type: 'static', steps: []});
    const forkIcon = (node: {root: HTMLElement}, isReadonly: boolean) =>
      node.root.querySelector(`dg-icon-fa[aria-label*="${sourceRunTooltip(isReadonly)}"]`);
    try {
      expect(!!forkIcon(basedOn, false)?.querySelector('.fa-code-branch'), true);
      expect(!!basedOn.root.querySelector('dg-icon-fa[aria-label="Remove"]'), true);
      expect(forkIcon(basedOn, false)?.nextElementSibling, null);
      expect(!!basedOnStep.root.querySelector(`dg-icon-fa[aria-label="${sourceRunTooltip(false, 'Baseline')}"]`), true);
      expect(basedOn.root.querySelector('.fa-lock'), null);
      for (const node of [asIs, asIsPipeline]) {
        expect(!!forkIcon(node, true)?.querySelector('.fa-code-branch'), true);
        expect(!!forkIcon(node, true)?.closest('.ml-auto'), true);
        const lock = node.root.querySelector('.fa-lock')?.closest('dg-icon-fa');
        expect(lock?.getAttribute('aria-label')?.includes(sourceRunLabel(true)), true);
        expect(lock?.getAttribute('aria-label')?.includes('Click for options'), false);
        expect(!!lock?.closest('.ml-auto'), false);
        expect(node.root.querySelectorAll('dg-icon-fa .fa-lock').length, 1);
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

  test('The source run menu offers comparison for editable nodes only', async () => {
    // the menu renders on the next tick; its first label is the header
    const menuItems = async (node: {root: HTMLElement}) => {
      node.root.querySelector('.fa-code-branch')!.closest<HTMLElement>('dg-icon-fa')!.click();
      await waitUntil(() => document.querySelector('.d4-menu-popup') !== null, 'source run menu not shown', 3000);
      const items = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
        .map((el) => el.textContent?.trim()).slice(1);
      document.querySelectorAll('.d4-menu-popup').forEach((el) => el.remove());
      return items;
    };
    const editable = mountNode({type: 'funccall', sourceRunId: 'run1'});
    const readonly = mountNode({type: 'funccall', isReadonly: true, sourceRunId: 'run2'}, {isReadonly: true});
    try {
      expectArray(await menuItems(editable), ['Open original run', 'Compare with original', 'Add to comparison...']);
      expectArray(await menuItems(readonly), ['Open original run']);
      expect(readonly.root.querySelector('.fa-lock')!.closest('dg-icon-fa')!.getAttribute('role'), 'img');
    } finally {
      document.querySelectorAll('.d4-menu-popup').forEach((el) => el.remove());
      editable.unmount();
      readonly.unmount();
    }
  });

  test('The selected node shows its actions with labels', async () => {
    const selected = mountNode({type: 'funccall'}, {isOutputOutdated: true});
    const other = mountNode({type: 'funccall'}, {isSelected: false, isOutputOutdated: true});
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
