import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';
import {initSearch, powerSearch, queriesSearch} from "./search/power-search";
import {widgetHost, getSettings, UserWidgetsSettings } from './utils';

let settings: UserWidgetsSettings;

export function welcomeView() {

  let searchStr = null;
  if (window.location.pathname == '/search' && window.location.search.startsWith('?')) {
    const params = new URLSearchParams(window.location.search.slice(1));
    searchStr = params.get('q') ?? '';
  }

  // let's not add this view if the platform was opened via a specific link
  if (searchStr == null && window.location.pathname != '/' && window.location.pathname != '/login.html')
    return;

  let input = ui.element('input', 'ui-input-editor') as HTMLInputElement;
  input.placeholder = 'Search everywhere. Try "aspirin" or "7JZK"';
  let inputHost = ui.div([
    ui.iconFA('search'),
    ui.div([
      input
    ], 'ui-input-root,ui-input-type-ahead')
  ], 'd4-search-bar');

  let searchHost = ui.block([], 'power-pack-search-host');
  let widgetsHost = ui.div([], 'power-pack-widgets-host');
  let viewHost = ui.div([widgetsHost, searchHost]);
  let view = grok.shell.newView(
    'Welcome',
    [inputHost, viewHost],
    'power-pack-welcome-view');

  let widgetFunctions = DG.Func.find({tags: ['dashboard'], returnType: 'widget'});

  getSettings().then((settings: UserWidgetsSettings) => {
    for (let f of widgetFunctions) {
      if (!settings[f.name] || settings[f.name].ignored) {
        const widgetHeader = ui.div();
        f.apply({'header': widgetHeader}).then(function (w: DG.Widget) {
          w.factory = f;
          widgetsHost.appendChild(widgetHost(w, widgetHeader));
        }).catch((e) => {
          console.error(`Unable to execute function ${f.name}`, e);
        });
      }
    }
  });

  initSearch();

  function doSearch(s: string) {
    input.value = s;
    let search = s !== '';
    widgetsHost.style.display = (search ? 'none' : '');
    searchHost.style.display = (search ? '' : 'none');
    if (search != null)
      powerSearch(s, searchHost);
    view.path = search ? `search?q=${s}` : 'search';
  }

  rxjs.fromEvent(input, 'input').pipe(debounceTime(500)).subscribe(_ => doSearch(input.value));

  if (searchStr != null)
    doSearch(searchStr);
}
