import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {entityIcon} from './workspace-tab';
import {
  querySearch, scriptsSearch,
  functionSearch, appSearch,
} from '../search/entity-search';


const TYPE_OPTIONS = ['All', 'Projects', 'Connections', 'Queries', 'Scripts', 'Functions', 'Apps',
  'Layouts', 'Notebooks', 'Models'] as const;

type EntityType = typeof TYPE_OPTIONS[number];

/** Interleaves arrays round-robin so every category is represented. */
function interleave(arrays: DG.Entity[][]): DG.Entity[] {
  const result: DG.Entity[] = [];
  const seen = new Set<string>();
  let i = 0;
  let added = true;
  while (added) {
    added = false;
    for (const arr of arrays) {
      if (i < arr.length) {
        const e = arr[i];
        if (!seen.has(e.id)) {
          seen.add(e.id);
          result.push(e);
        }
        added = true;
      }
    }
    i++;
  }
  return result;
}

const MAX_RESULTS = 200;

async function searchEntities(query: string, type: EntityType): Promise<DG.Entity[]> {
  const dapiList = (src: any) => src.filter(query).by(MAX_RESULTS).list() as Promise<DG.Entity[]>;
  const clientList = (fn: (s: string) => Promise<DG.Entity[]>) =>
    fn(query).then((r) => r.slice(0, MAX_RESULTS));

  if (type === 'All') {
    const groups = await Promise.all([
      clientList(appSearch).catch(() => []),
      dapiList(grok.dapi.projects).catch(() => []),
      dapiList(grok.dapi.connections).catch(() => []),
      clientList(querySearch).catch(() => []),
      clientList(scriptsSearch).catch(() => []),
      clientList(functionSearch).catch(() => []),
      dapiList(grok.dapi.layouts).catch(() => []),
      dapiList(grok.dapi.notebooks).catch(() => []),
      dapiList(grok.dapi.models).catch(() => []),
    ]);
    return interleave(groups);
  }

  switch (type) {
  case 'Projects': return dapiList(grok.dapi.projects);
  case 'Connections': return dapiList(grok.dapi.connections);
  case 'Queries': return clientList(querySearch);
  case 'Scripts': return clientList(scriptsSearch);
  case 'Functions': return clientList(functionSearch);
  case 'Apps': return clientList(appSearch);
  case 'Layouts': return dapiList(grok.dapi.layouts);
  case 'Notebooks': return dapiList(grok.dapi.notebooks);
  case 'Models': return dapiList(grok.dapi.models);
  }
}


export function showManageFavoritesDialog(
  group: DG.Group,
  currentEntities: DG.Entity[],
  onChange: () => void,
): void {
  const pinnedIds = new Set(currentEntities.map((e) => e.id));

  // Pinned section
  const pinnedCountSpan = ui.span([`${currentEntities.length}`]);
  const pinnedHeader = ui.div(
    [ui.span(['Pinned (', pinnedCountSpan, '):'])],
    'pp-manage-fav-section-header',
  );
  const pinnedListDiv = ui.divV([], 'pp-manage-fav-list');
  const pinnedEmptyHint = ui.divText('No pinned items', 'pp-manage-fav-empty');

  function refreshPinnedList(entities: DG.Entity[]): void {
    pinnedListDiv.innerHTML = '';
    pinnedCountSpan.textContent = `${entities.length}`;
    if (entities.length === 0) {
      pinnedListDiv.appendChild(pinnedEmptyHint);
      return;
    }
    for (const entity of entities) {
      const icon = entityIcon(entity);
      icon.classList.add('pp-manage-fav-row-icon');
      const name = ui.div([entity.friendlyName], 'pp-manage-fav-row-name');
      const removeBtn = ui.iconFA('times', async () => {
        await DG.Favorites.remove(entity, group);
        pinnedIds.delete(entity.id);
        currentEntities = currentEntities.filter((e) => e.id !== entity.id);
        refreshPinnedList(currentEntities);
        onChange();
      });
      removeBtn.classList.add('pp-manage-fav-row-action');
      ui.tooltip.bind(removeBtn, 'Remove');
      pinnedListDiv.appendChild(ui.divH([icon, name, removeBtn], 'pp-manage-fav-row'));
    }
  }

  refreshPinnedList(currentEntities);

  // Search section
  const searchInput = ui.input.search('', {onValueChanged: () => triggerSearch()});
  searchInput.root.classList.add('pp-manage-fav-search-input');
  const typeInput = ui.input.choice('', {items: TYPE_OPTIONS as unknown as string[], value: 'All',
    onValueChanged: () => triggerSearch()});
  typeInput.root.classList.add('pp-manage-fav-type-input');
  const searchResultsDiv = ui.divV([], 'pp-manage-fav-list pp-manage-fav-search-results');
  const searchHint = ui.divText('Type to search...', 'pp-manage-fav-search-hint');
  searchResultsDiv.appendChild(searchHint);

  let searchTimer: ReturnType<typeof setTimeout> | null = null;
  let searchGeneration = 0;

  function triggerSearch(): void {
    if (searchTimer)
      clearTimeout(searchTimer);
    searchTimer = setTimeout(() => doSearch(), 300);
  }

  async function addEntity(entity: DG.Entity, row: HTMLElement): Promise<void> {
    await DG.Favorites.add(entity, group);
    pinnedIds.add(entity.id);
    currentEntities.push(entity);
    refreshPinnedList(currentEntities);
    row.remove();
    if (searchResultsDiv.querySelectorAll('.pp-manage-fav-row').length === 0)
      searchResultsDiv.appendChild(ui.divText('No results', 'pp-manage-fav-search-hint'));
    onChange();
  }

  async function doSearch(): Promise<void> {
    const query = (searchInput.value as string ?? '').trim();
    if (query.length < 2) {
      searchResultsDiv.innerHTML = '';
      searchResultsDiv.appendChild(searchHint);
      return;
    }

    const gen = ++searchGeneration;
    searchResultsDiv.innerHTML = '';
    searchResultsDiv.appendChild(ui.loader());

    try {
      const type = typeInput.value as EntityType;
      let results = await searchEntities(query, type);
      if (gen !== searchGeneration)
        return;
      results = results.filter((e) => !pinnedIds.has(e.id));
      searchResultsDiv.innerHTML = '';

      if (results.length === 0) {
        searchResultsDiv.appendChild(ui.divText('No results', 'pp-manage-fav-search-hint'));
        return;
      }

      for (const entity of results) {
        const icon = entityIcon(entity);
        icon.classList.add('pp-manage-fav-row-icon');
        const name = ui.div([entity.friendlyName], 'pp-manage-fav-row-name');
        const addBtn = ui.iconFA('plus', () => {});
        addBtn.classList.add('pp-manage-fav-row-action', 'pp-manage-fav-row-add');
        ui.tooltip.bind(addBtn, 'Pin');
        const row = ui.divH([icon, name, addBtn], 'pp-manage-fav-row pp-manage-fav-row-clickable');
        row.addEventListener('click', () => addEntity(entity, row));
        searchResultsDiv.appendChild(row);
      }
    }
    catch (e) {
      if (gen !== searchGeneration)
        return;
      searchResultsDiv.innerHTML = '';
      searchResultsDiv.appendChild(ui.divText('Search failed', 'pp-manage-fav-search-hint'));
      console.error('Manage favorites search failed', e);
    }
  }

  const searchBar = ui.divH([searchInput.root, typeInput.root], 'pp-manage-fav-search-bar');
  const addLabel = ui.div(['Add:'], 'pp-manage-fav-section-header');

  const content = ui.divV([
    pinnedHeader, pinnedListDiv,
    addLabel, searchBar, searchResultsDiv,
  ], 'pp-manage-fav-content');

  const dlg = ui.dialog({title: `Manage Favorites: ${group.friendlyName || group.name}`});
  dlg.add(content);
  dlg.addButton('Close', () => dlg.close());
  // Hide default OK/Cancel
  const footer = dlg.root.querySelector('.d4-dialog-footer');
  if (footer) {
    const buttons = footer.querySelectorAll('.ui-btn');
    for (const btn of Array.from(buttons))
      if (btn.textContent === 'OK' || btn.textContent === 'CANCEL')
        (btn as HTMLElement).style.display = 'none';
  }

  dlg.show({resizable: true, width: 480, height: 520});
}
