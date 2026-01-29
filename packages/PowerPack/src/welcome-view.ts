/* eslint-disable max-len */
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
  input.classList.add('power-search-search-everywhere-input'); // otherwise HTML complains if I do this above :)
  input.placeholder = 'Search everywhere. Try "aspirin" or "7JZK"';
  const inputContainer = ui.div([
    input,
  ], 'ui-input-root,ui-input-type-ahead');
  const inputHost = ui.div([
    ui.iconFA('search'),
    inputContainer,
  ], 'd4-search-bar');
  const helpDiv = ui.divText('Use Ctrl+Space to see suggestions', {style: {
    marginLeft: '12px', fontSize: '10px', color: 'var(--text-color-light)',
    fontStyle: 'italic',
    marginBottom: '10px', flexGrow: '0',
  }, classes: 'power-search-help-text-container'});

  // input.addEventListener('input', () => {
  //   if (!input.value)
  //     helpDiv.style.visibility = 'hidden';
  //   else
  //     helpDiv.style.visibility = 'visible';
  // });

  // helpDiv.style.display = input.value ? 'flex' : 'none';
  inputHost.style.marginBottom = '0px';
  suggestionMenuKeyNavigation(inputContainer);

  const searchHost = ui.block([], 'power-pack-search-host');
  const widgetsHost = ui.div([], 'power-pack-widgets-host');
  const widgetsPanel = ui.div([widgetsHost]);
  const viewHost = ui.div([widgetsPanel, searchHost]);
  const view = DG.View.create();
  view.root.appendChild(inputHost);
  view.root.appendChild(helpDiv);
  view.root.appendChild(viewHost);
  view.root.classList.add('power-pack-welcome-view');

  const widgetFunctions = DG.Func.find({meta: {role: DG.FUNC_TYPES.DASHBOARD}, returnType: 'widget'});
  widgetFunctions.sort((a, b) => (a.options['order'] ?? 100) - (b.options['order'] ?? 100));
  const widgetHosts: {[index: string]: HTMLElement} = {};
  const settings: UserWidgetsSettings = getSettings();

  function refresh() {
    grok.dapi.groups.find(DG.User.current().group.id).then((userGroup: DG.Group) => {
      while (widgetsHost.firstChild)
        widgetsHost.removeChild(widgetsHost.firstChild);

      for (const f of widgetFunctions) {
        const canView: string[] = f.options['canView']?.split(',') ?? [];
        if (canView.length === 0 || (userGroup.memberships.some((g) => canView.includes(g.friendlyName)) ||
            userGroup.adminMemberships.some((g) => canView.includes(g.friendlyName)))) {
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
        },
      })),
    );
  }

  widgetsPanel.appendChild(ui.link('Customize widgets...', () => customizeWidgets()));

  function doSearch(s: string) {
    input.value = s;
    const search = s !== '';
    widgetsPanel.style.display = (search ? 'none' : '');
    searchHost.style.display = (search ? '' : 'none');
    if (search != null)
      powerSearch(s, searchHost, input);
    view.path = search ? `search?q=${encodeURIComponent(s)}` : 'search';
  }

  rxjs.fromEvent(input, 'input').pipe(debounceTime(500)).subscribe((_) => doSearch(input.value));

  if (searchStr != null)
    doSearch(searchStr);
  return view;
}

function suggestionMenuKeyNavigation(inputContainer: HTMLElement) {
  inputContainer.addEventListener('keydown', (e) => {
    if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp' && e.key !== 'Enter')
      return;
    const currentlySelected: HTMLElement | null = inputContainer.querySelector('.d4-menu-item-hover');
    const allItems: HTMLElement[] = Array.from(inputContainer.querySelectorAll('.d4-menu-item') ?? []);
    if (!allItems || allItems.length === 0)
      return;
    allItems.sort((a, b) => a.offsetTop - b.offsetTop); // sort by vertical position

    let currentIndex = currentlySelected ? allItems.indexOf(currentlySelected) : -1;
    if (e.key === 'ArrowDown')
      currentIndex = (currentIndex + 1) % allItems.length;
    else if (e.key === 'ArrowUp')
      currentIndex = (currentIndex - 1 + allItems.length) % allItems.length;
    else if (e.key === 'Enter' && currentlySelected) {
      currentlySelected.click();
      e.preventDefault();
      return;
    }
    currentlySelected?.classList?.remove('d4-menu-item-hover');
    allItems[currentIndex]?.classList?.add('d4-menu-item-hover');
  });
}
