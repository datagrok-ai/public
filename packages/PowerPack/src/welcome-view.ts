import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

function card(w: DG.Widget): HTMLElement {
  let host = ui.divV([], 'power-pack-widget-host');

  let header = ui.div([
    ui.divText(w.props.caption ?? '', 'd4-dialog-title'),
    ui.icons.settings(() => { grok.shell.o = w}, 'Edit settings'),
    ui.icons.close(() => host.remove(), 'Remove'),
  ], 'd4-dialog-header');

  ui.tools.setHoverVisibility(host, Array.from(host.querySelectorAll('i')));

  host.appendChild(header);
  host.appendChild(ui.div([w.root], 'power-pack-widget-content'));
  return host;
}

export function welcomeView() {
  let view = grok.shell.newView('Welcome');

  let widgetFunctions = DG.Func.find({returnType: 'widget'});

  for (let f of widgetFunctions)
    f.apply().then((w: DG.Widget) => view.root.appendChild(card(w)));
}