import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';
import {powerSearch} from './search/power-search';
import {widgetHost, getSettings, UserWidgetsSettings} from './utils';

export function welcomeView(): DG.View | undefined {
  let searchStr = null;
  if (window.location.pathname == '/search' && window.location.search.startsWith('?')) {
    const params = new URLSearchParams(window.location.search.slice(1));
    searchStr = params.get('q') ?? '';
  }

  const input = ui.element('input', 'ui-input-editor') as HTMLInputElement;
  input.placeholder = 'Search everywhere. Try "aspirin" or "7JZK"';
  const inputHost = ui.div([
    ui.iconFA('search'),
    ui.div([
      input,
    ], 'ui-input-root,ui-input-type-ahead'),
  ], 'd4-search-bar');

  const searchHost = ui.block([], 'power-pack-search-host');
  const widgetsHost = ui.div([], 'power-pack-widgets-host');
  const viewHost = ui.div([widgetsHost, searchHost]);
  const view = DG.View.create();
  view.root.appendChild(inputHost);
  view.root.appendChild(viewHost);
  view.root.classList.add('power-pack-welcome-view');

  const widgetFunctions = DG.Func.find({tags: ['dashboard'], returnType: 'widget'});
  const settings: UserWidgetsSettings = getSettings();
  for (const f of widgetFunctions) {
    if (!settings[f.name] || settings[f.name].ignored) {
      const widgetHeader = ui.div();
      f.apply({'header': widgetHeader}).then(function(w: DG.Widget) {
        if (!w)
          return;
        w.factory = f;
        widgetsHost.appendChild(widgetHost(w, widgetHeader));
      }).catch((e) => {
        console.error(`Unable to execute function ${f.name}`, e);
      });
    }
  }

  function doSearch(s: string) {
    input.value = s;
    const search = s !== '';
    widgetsHost.style.display = (search ? 'none' : '');
    searchHost.style.display = (search ? '' : 'none');
    if (search != null)
      powerSearch(s, searchHost);
    view.path = search ? `search?q=${s}` : 'search';
  }

  rxjs.fromEvent(input, 'input').pipe(debounceTime(500)).subscribe((_) => doSearch(input.value));

  if (searchStr != null)
    doSearch(searchStr);
  return view;
}
