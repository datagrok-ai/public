/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const WIDGETS_STORAGE = 'widgets';

export function card(w: DG.Widget): HTMLElement {
  let host = ui.box(null, 'power-pack-widget-host');

  function remove(): void {
    host.remove();
    if (w.factory?.name) {
      let settings = { ignored: true };
      grok.dapi.userDataStorage
        .postValue(WIDGETS_STORAGE, w.factory.name, JSON.stringify(settings))
        .then((_) => grok.shell.info('To control widget visibility, go to Tools | Widgets'));
    }
  }

  let header = ui.div([
    ui.divText(w.props.caption ?? '', 'd4-dialog-title'),
    ui.icons.settings(() => { grok.shell.o = w}, 'Edit settings'),
    ui.icons.close(remove, 'Remove'),
  ], 'd4-dialog-header');

  host.appendChild(header);
  host.appendChild(ui.box(w.root, 'power-pack-widget-content'));
  ui.tools.setHoverVisibility(host, Array.from(host.querySelectorAll('i')));

  return host;
}