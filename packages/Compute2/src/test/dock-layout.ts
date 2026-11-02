import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, waitUntil, delay, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import {awaitWebComponents, closeView} from './utils';

const HOST_WIDTH = 1000;

// A background view is mounted while hidden, so its dock gets a size only when shown
function hiddenHost() {
  const host = document.createElement('div');
  host.style.cssText = `position:fixed;left:0;top:0;width:${HOST_WIDTH}px;height:600px;display:none`;
  document.body.append(host);
  return host;
}

const widthRatio = (el: Element | null) =>
  el ? el.getBoundingClientRect().width / HOST_WIDTH : NaN;

const awaitNarrowForm = (form: () => Element | null) => waitUntil(() => Math.abs(widthRatio(form()) - 0.2) < 0.03,
  'form should take ~20% of the width', 5000);

const titles = (dock: Element) => dock.shadowRoot?.textContent ?? '';

async function mountDock(host: HTMLElement) {
  const dock = document.createElement('dock-spawn-ts');
  dock.style.cssText = 'display:block;width:100%;height:100%';
  const initialized = new Promise((resolve) =>
    dock.addEventListener('manager-init-finished', resolve, {once: true}));
  host.append(dock);
  await initialized;
  const panel = (title: string) => {
    const el = document.createElement('div');
    el.setAttribute('dock-spawn-title', title);
    return el;
  };
  return {dock, panel};
}

category('Dock: layout of hidden mounts', () => {
  test('Docking ratio is kept when the dock is shown later', async () => {
    await awaitWebComponents();
    const host = hiddenHost();
    try {
      const {dock, panel} = await mountDock(host);
      const form = panel('Inputs');
      form.setAttribute('dock-spawn-dock-type', 'left');
      form.setAttribute('dock-spawn-dock-ratio', '0.2');
      dock.append(form);
      // one task, so the dock sees the form while still hidden
      await delay(0);
      host.style.display = 'block';
      await awaitNarrowForm(() => form);
    } finally {
      host.remove();
    }
  });

  test('Docking ratio is kept when a panel is replaced while the dock is hidden', async () => {
    await awaitWebComponents();
    const host = hiddenHost();
    host.style.display = 'block';
    try {
      const {dock, panel} = await mountDock(host);
      const form = panel('Steps');
      form.setAttribute('dock-spawn-dock-type', 'left');
      form.setAttribute('dock-spawn-dock-ratio', '0.2');
      const review = panel('Review');
      dock.append(review, form);
      await awaitNarrowForm(() => form);
      // a view in the background re-renders its main panel, as a workflow does when its run loads
      host.style.display = 'none';
      await delay(0);
      review.remove();
      dock.append(panel('Review'));
      await delay(0);
      host.style.display = 'block';
      await awaitNarrowForm(() => form);
    } finally {
      host.remove();
    }
  });

  test('A panel removed while the dock is hidden is gone once it is shown', async () => {
    await awaitWebComponents();
    const host = hiddenHost();
    host.style.display = 'block';
    try {
      const {dock, panel} = await mountDock(host);
      const form = panel('Steps');
      form.setAttribute('dock-spawn-dock-type', 'left');
      form.setAttribute('dock-spawn-dock-ratio', '0.2');
      const extra = panel('Extra panel');
      dock.append(panel('Review'), extra, form);
      await awaitCheck(() => titles(dock).includes('Extra panel'), 'the extra panel should dock', 5000);
      host.style.display = 'none';
      await delay(0);
      extra.remove();
      await delay(0);
      host.style.display = 'block';
      await awaitCheck(() => !titles(dock).includes('Extra panel'), 'the removed panel should close', 5000);
      await awaitNarrowForm(() => form);
    } finally {
      host.remove();
    }
  });

  test('A panel added and removed while the dock is hidden never docks', async () => {
    await awaitWebComponents();
    const host = hiddenHost();
    try {
      const {dock, panel} = await mountDock(host);
      const review = panel('Review');
      const extra = panel('Extra panel');
      dock.append(review, extra);
      await delay(0);
      extra.remove();
      await delay(0);
      host.style.display = 'block';
      await awaitCheck(() => titles(dock).includes('Review'), 'the review panel should dock', 5000);
      expect(titles(dock).includes('Extra panel'), false);
    } finally {
      host.remove();
    }
  });

  test('RFV mounted in the background keeps the narrow form', async () => {
    await awaitWebComponents();
    const call = DG.Func.byName('Compute2:GraphicsSvgOutputTest').prepare({a: 1});
    const view = await grok.functions.call('Compute2:RichFunctionViewEditor', {call}) as DG.ViewBase;
    const host = hiddenHost();
    try {
      view.root.style.height = '100%';
      host.append(view.root);
      const inputs = () => view.root.querySelector('[dock-spawn-title="Inputs"]');
      await waitUntil(() => inputs()?.querySelector('dg-input-form') != null, 'RFV form not rendered', 15000);
      host.style.display = 'block';
      await awaitNarrowForm(inputs);
    } finally {
      host.remove();
      closeView(view);
    }
  });
});
