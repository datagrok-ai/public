import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function welcomeView() {
  let view = grok.shell.newView('Welcome');

  let widgetFunctions = DG.Func.find({returnType: 'widget'});

  for (let f of widgetFunctions)
    f.apply().then((w: DG.Widget) => view.root.appendChild(w.root));
}