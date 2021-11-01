/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
//import { ODEview } from './views/ode_view';

export const _package = new DG.Package();

//name: SimPKPD
//tags: app
export async function sim() {
  grok.shell.info("Hello world");

  let result = await grok.functions.call(
    "Simpkpd:rxodeCommandReal", {
    "inputSD": 1
  });

  const view = grok.shell.newView('Custom View', [ui.image(`data:image/png;base64,${result}`, 500, 500)]);
}
