import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';
import {powerSearch} from './search/power-search';
import {getSettings, saveSettings, UserWidgetsSettings, widgetHostFromFunc} from './utils';

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
  const widgetsPanel = ui.div([widgetsHost]);
  const viewHost = ui.div([widgetsPanel, searchHost]);
  const view = DG.View.create();
  view.root.appendChild(inputHost);
  view.root.appendChild(viewHost);
  view.root.classList.add('power-pack-welcome-view');

  const widgetFunctions = DG.Func.find({tags: ['dashboard'], returnType: 'widget'});
  const widgetHosts: {[index: string]: HTMLElement} = {};
  const settings: UserWidgetsSettings = getSettings();

  function refresh() {
    grok.dapi.groups.find(DG.User.current().group.id).then((userGroup: DG.Group) => {
      while (widgetsHost.firstChild)
        widgetsHost.removeChild(widgetsHost.firstChild);

      for (const f of widgetFunctions) {
        const canView: string[] = f.options['canView']?.split(',') ?? [];
        if (canView.length === 0 || (userGroup.memberships.some((g) => canView.includes(g.friendlyName))
            || userGroup.adminMemberships.some((g) => canView.includes(g.friendlyName)))) {
          if (!settings[f.name] || !settings[f.name].ignored)
            widgetsHost.appendChild(widgetHosts[f.name] ??= widgetHostFromFunc(f));
        }
      }
    });
  }

  refresh();

  function customizeWidgets() {
    grok.shell.windows.context.visible = true;
    const existingNames = Object.keys(settings).filter((name) => DG.Func.byName(name));

    grok.shell.o = ui.form(
        existingNames.map((name) => ui.input.bool(DG.Func.byName(name).friendlyName, {
          value: !settings[name].ignored,
          onValueChanged: (value, input) => {
            settings[name].ignored = !value;
            refresh();
            saveSettings();
          }
        }))
    );
  }

  widgetsPanel.appendChild(ui.link('Customize widgets...', () => customizeWidgets()));

  function doSearch(s: string) {
    input.value = s;
    const search = s !== '';
    widgetsPanel.style.display = (search ? 'none' : '');
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
