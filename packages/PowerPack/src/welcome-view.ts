/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';
import {powerSearch, showSearchSuggestions} from './search/power-search';
import {getSettings, saveSettings, UserWidgetsSettings, widgetHostFromFunc} from './utils';
import {setWorkspacePreviewHost, setWorkspaceCustomizeEl} from './spotlight/preview-host';
import {getCurrentUserGroup} from './spotlight/group-favorites';

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
  ], 'd4-search-bar,power-pack-welcome-search');
  const helpDiv = ui.divText('Use Ctrl+Space to see suggestions', 'power-search-help-text-container');
  inputHost.appendChild(helpDiv);

  const updateHelp = () => {helpDiv.style.visibility = input.value ? 'visible' : 'hidden';};
  updateHelp();
  input.addEventListener('input', updateHelp);
  input.addEventListener('keydown', (e) => {
    if (e.ctrlKey && e.code === 'Space') {
      e.preventDefault();
      showSearchSuggestions(input);
    }
  });
  suggestionMenuKeyNavigation(inputContainer);

  const searchHost = ui.block([], 'power-pack-search-host');
  searchHost.setAttribute('aria-busy', 'false');
  const widgetsHost = ui.div([], 'power-pack-widgets-host');
  const previewHost = ui.div([], 'power-pack-workspace-preview');
  previewHost.style.display = 'none';
  const widgetsPanel = ui.div([widgetsHost, previewHost], 'power-pack-widgets-panel');
  const viewHost = ui.div([widgetsPanel, searchHost], 'power-pack-view-host');
  setWorkspacePreviewHost(previewHost);
  const view = DG.View.create();
  view.root.appendChild(viewHost);
  // The search takes the otherwise empty top row: the view's ribbon without tabs, or the tab strip while Home is its
  // only tab. The mode is known only after the settings load.
  let place: 'ribbon' | 'strip' | 'page' | null = null;
  const tabStrip = () => view.root.closest('.grok-tab-host')?.querySelector<HTMLElement>(':scope > .tab-handle-list-container');
  const stripRoom = (strip: HTMLElement) => {
    let used = 0;
    for (const c of Array.from(strip.children))
      if (c !== inputHost) used += (c as HTMLElement).offsetWidth;
    return strip.clientWidth - used;
  };
  const placeSearch = () => {
    const strip = tabStrip();
    const newPlace = grok.shell.windows.simpleMode ? 'ribbon' :
      strip && strip.querySelectorAll('.tab-handle:not(.tab-handle-preview, .tab-handle-temp-tab)').length === 1 &&
      stripRoom(strip) >= 300 ? 'strip' : 'page';
    if (newPlace === place)
      return;
    place = newPlace;
    view.setRibbonPanels(place === 'ribbon' ? [[inputHost]] : []);
    if (place === 'strip')
      strip!.append(inputHost);
    else if (place === 'page')
      view.root.prepend(inputHost);
  };
  placeSearch();
  const resized = new rxjs.Subject<void>();
  const resizeObserver = new ResizeObserver(() => resized.next());
  resizeObserver.observe(view.root);
  view.subs.push(
    rxjs.merge(grok.events.onEvent('grok-panels-changed'), grok.events.onViewAdded, grok.events.onViewRemoved, resized).pipe(debounceTime(50)).subscribe(placeSearch),
    new rxjs.Subscription(() => resizeObserver.disconnect()));
  view.root.classList.add('power-pack-welcome-view');

  const widgetFunctions = DG.Func.find({meta: {role: DG.FUNC_TYPES.DASHBOARD}, returnType: 'widget'});
  widgetFunctions.sort((a, b) => (a.options['order'] ?? 100) - (b.options['order'] ?? 100));
  const widgetHosts: {[index: string]: HTMLElement} = {};
  const settings: UserWidgetsSettings = getSettings();

  function refresh() {
    getCurrentUserGroup().then((userGroup) => {
      const userGroups = [...userGroup?.memberships ?? [], ...userGroup?.adminMemberships ?? []];

      while (widgetsHost.firstChild)
        widgetsHost.removeChild(widgetsHost.firstChild);

      for (const f of widgetFunctions) {
        const canView: string[] = f.options['canView']?.split(',') ?? [];
        if (canView.length === 0 || userGroups.some((g) => canView.includes(g.friendlyName))) {
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

  const customizeLink = ui.link('Customize widgets...', () => customizeWidgets());
  widgetsPanel.appendChild(customizeLink);
  setWorkspaceCustomizeEl(customizeLink);

  function doSearch(s: string) {
    input.value = s;
    updateHelp();
    const search = s !== '';
    widgetsPanel.style.display = (search ? 'none' : '');
    searchHost.style.display = (search ? '' : 'none');
    if (search != null)
      powerSearch(s, searchHost, input);
    view.path = search ? `search?q=${encodeURIComponent(s)}` : 'search';
  }

  input.addEventListener('input', () => searchHost.setAttribute('aria-busy', 'true'));
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
