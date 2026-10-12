import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
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

const awaitNarrowForm = (form: () => Element | null) => awaitCheck(() => Math.abs(widthRatio(form()) - 0.2) < 0.03,
  'form should take ~20% of the width', 5000);

category('Dock: layout of hidden mounts', () => {
  test('Docking ratio is kept when the dock is shown later', async () => {
    await awaitWebComponents();
    const host = hiddenHost();
    try {
      const dock = document.createElement('dock-spawn-ts');
      dock.style.cssText = 'display:block;width:100%;height:100%';
      const initialized = new Promise((resolve) =>
        dock.addEventListener('manager-init-finished', resolve, {once: true}));
      host.append(dock);
      await initialized;
      const form = document.createElement('div');
      form.setAttribute('dock-spawn-title', 'Inputs');
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

  test('RFV mounted in the background keeps the narrow form', async () => {
    await awaitWebComponents();
    const call = DG.Func.byName('Compute2:GraphicsSvgOutputTest').prepare({a: 1});
    const view = await grok.functions.call('Compute2:RichFunctionViewEditor', {call}) as DG.ViewBase;
    const host = hiddenHost();
    try {
      view.root.style.height = '100%';
      host.append(view.root);
      const inputs = () => view.root.querySelector('[dock-spawn-title="Inputs"]');
      await awaitCheck(() => inputs()?.querySelector('dg-input-form') != null, 'RFV form not rendered', 15000);
      host.style.display = 'block';
      await awaitNarrowForm(inputs);
    } finally {
      host.remove();
      closeView(view);
    }
  });
});
