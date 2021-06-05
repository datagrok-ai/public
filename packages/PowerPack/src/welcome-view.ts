import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';

interface UserWidgetSettings {
  factoryName?: string;
  caption?: string;
  ignored?: boolean;
}

interface UserWidgetsSettings {
  [index: string]: UserWidgetSettings;
}

const WIDGETS_STORAGE = 'widgets';
let settings: UserWidgetsSettings;

function card(w: DG.Widget): HTMLElement {
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

export function welcomeView() {
  let input = ui.element('input', 'ui-input-editor') as HTMLInputElement;
  input.placeholder = 'Search everywhere. Try "aspirin" or "7JZK"';
  let inputHost = ui.div([
    ui.iconFA('search'),
    ui.div([
      input
    ], 'ui-input-root,ui-input-type-ahead')
  ], 'd4-search-bar');

  let searchHost = ui.divV([], 'power-pack-search-host');
  let widgetsHost = ui.div([], 'power-pack-widgets-host');
  let viewHost = ui.div([widgetsHost, searchHost]);
  grok.shell.newView('Welcome', [inputHost, viewHost], 'power-pack-welcome-view');

  let widgetFunctions = DG.Func.find({returnType: 'widget'});
  let searchFunctions = DG.Func.find({tags: ['search'], returnType: 'list'});
  let searchWidgetFunctions = DG.Func.find({tags: ['search'], returnType: 'widget'});

  grok.dapi.userDataStorage.get(WIDGETS_STORAGE).then((settings) => {
    for (let f of widgetFunctions) {
      if (!settings[f.name] || settings[f.name].ignored)
        f.apply().then(function (w: DG.Widget) {
          w.factory = f;
          widgetsHost.appendChild(card(w));
        });
    }
  });

  function doSearch(s: string) {
    ui.empty(searchHost);

    if (DG.View.ALL_VIEW_TYPES.includes(s)) {
      searchHost.appendChild(DG.View.createByType(s).root);
      return;
    }

    for (let sf of searchFunctions)
      sf.apply({s: input.value}).then((results: any[]) => {
        if (results.length > 0) {
          searchHost.appendChild(ui.divV([
            ui.h3(sf.description ?? sf.name),
            ui.list(results)
          ]));
        }
    });

    grok.functions
      .eval(s)
      .then((result) => searchHost.appendChild(ui.span([s + ' = ', result], {style: {'font-size' : '20px'}})))
      .catch(() => {})

    for (let sf of searchWidgetFunctions)
      sf.apply({s: input.value})
        .then((result: DG.Widget) => {
          if (result)
            searchHost.appendChild(card(result));
        });
  }

  rxjs.fromEvent(input, 'input').pipe(debounceTime(300)).subscribe(_ => {
    let search = input.value !== '';
    widgetsHost.style.display = (search ? 'none' : '');
    searchHost.style.display = (search ? '' : 'none');
    doSearch(input.value);
  });
}