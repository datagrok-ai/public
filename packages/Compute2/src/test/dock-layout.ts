import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, delay} from '@datagrok-libraries/test/src/test';
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

category('Dock: layout of hidden mounts', () => {
  test('Docking ratio is kept when the dock is shown later', async () => {
    await awaitWebComponents();
    const host = hiddenHost();
    try {
      const dock = document.createElement('dock-spawn-ts');
      dock.style.cssText = 'display:block;width:100%;height:100%';
      host.append(dock);
      await delay(100);
      const form = document.createElement('div');
      form.setAttribute('dock-spawn-title', 'Inputs');
      form.setAttribute('dock-spawn-dock-type', 'left');
      form.setAttribute('dock-spawn-dock-ratio', '0.2');
      dock.append(form);
      await delay(300);
      host.style.display = 'block';
      await delay(500);
      const ratio = widthRatio(form);
      expect(Math.abs(ratio - 0.2) < 0.03, true, `form should take ~20% of the width, got ${ratio.toFixed(2)}`);
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
      await delay(2000);
      host.style.display = 'block';
      await delay(1000);
      const ratio = widthRatio(view.root.querySelector('[dock-spawn-title="Inputs"]'));
      expect(Math.abs(ratio - 0.2) < 0.03, true, `form should take ~20% of the width, got ${ratio.toFixed(2)}`);
    } finally {
      host.remove();
      closeView(view);
    }
  });
});
