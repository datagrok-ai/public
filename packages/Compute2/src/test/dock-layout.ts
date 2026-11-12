import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, waitUntil} from '@datagrok-libraries/test/src/test';
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

category('Dock: layout of hidden mounts', () => {
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
