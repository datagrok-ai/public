import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {queries} from '../package-api';
import {UaFilter} from '../filter';
import {emptyState, formatGridTimes, formatTime, onRowContextMenu, rowsTable} from '../utils';
import {ErrorsView} from './errors';
import {debounceTime} from 'rxjs/operators';

const STALE_MS = 30000;

export class ClicksView extends UaView {
  expanded: {[key: string]: boolean} = {f: true, l: true};
  tabControl?: DG.TabControl;
  followedView: string = '';
  private refreshers: (() => void)[] = [];
  private clicksDf?: DG.DataFrame;

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Clicks';
  }

  async initViewers(path?: string): Promise<void> {
    this.root.className = 'grok-view ui-box';
    const tabs: {[key: string]: (() => HTMLElement)} = {
      'Click Analysis': () => this.createFilteredElement((filter) => this.getClickAnalysisTab(filter)),
      'Clicks': () => this.createFilteredElement((filter) => this.getClicksTab(filter)),
      'Followed by Error': () => this.createFilteredElement((filter) => this.getFollowedByErrorTab(filter)),
    };
    this.tabControl = ui.tabControl(tabs);
    this.root.appendChild(this.tabControl.root);
    const refresh = () => {
      for (const r of this.refreshers)
        r();
    };
    this.tabControl.onTabChanged.subscribe(() => {
      grok.shell.o = null;
      refresh();
    });
    this.uaToolbox.viewHandler.view.tabs.onTabChanged.subscribe(() => refresh());
  }

  async exportFiles(): Promise<DG.FileInfo[]> {
    return UaView.csvFiles({'clicks': this.clicksDf});
  }

  /** A sub-tab built for the applied filter, rebuilt when it is shown after the filter changed or
   * {@link STALE_MS} after it was built. */
  createFilteredElement(build: (filter: UaFilter) => Promise<HTMLElement>): HTMLElement {
    let shown = this.uaToolbox.filterStream.value;
    let current = this.createWaitElement(shown, build);
    let builtAt = Date.now();
    const refresh = () => {
      const filter = this.uaToolbox.filterStream.value;
      if (!current.isConnected || (filter === shown && Date.now() - builtAt < STALE_MS))
        return;
      shown = filter;
      builtAt = Date.now();
      const next = this.createWaitElement(filter, build);
      current.replaceWith(next);
      current = next;
    };
    this.refreshers.push(refresh);
    this.uaToolbox.filterStream.subscribe(() => refresh());
    return current;
  }

  createWaitElement(filter: UaFilter, build: (filter: UaFilter) => Promise<HTMLElement>): HTMLElement {
    const elem: HTMLElement = ui.wait(async () => build(filter));
    elem.style.height = '100%';
    return elem;
  }

  async getClicksTab(filter: UaFilter): Promise<HTMLElement> {
    const table = await queries.clicks(filter.date!, filter.groups);
    table.name = 'Clicks';
    const grid = DG.Viewer.grid(table, {showRowHeader: false, allowRowSelection: false, allowBlockSelection: false});
    grid.col('element')!.width = 400;
    formatGridTimes(grid);
    onRowContextMenu(grid, (menu, i) => {
      const action = table.get('action id', i);
      if (action)
        menu.item('Timeline', () => this.openTimeline('action', action));
    });
    table.onCurrentRowChanged.subscribe(() => {
      const i = table.currentRowIdx;
      if (i < 0)
        return;
      const action = table.get('action id', i);
      const acc = DG.Accordion.create();
      acc.addPane(table.get('type', i), () => ui.divV([
        ui.tableFromMap({
          'Time': formatTime(table.get('time', i)),
          'User': table.get('user', i) ?? '',
          'Type': table.get('type', i),
          'Element': table.get('element', i),
          ...(action ? {'Action id': action} : {}),
        }),
        action ? ui.buttonsInput([ui.button('Timeline', () => this.openTimeline('action', action))]) :
          ui.divText('No action id'),
      ]), true);
      grok.shell.o = acc.root;
    });
    return grid.root;
  }

  /** Clicks per element and how many of them an error followed within 5 s (same action id). */
  async getFollowedByErrorTab(filter: UaFilter): Promise<HTMLElement> {
    const view = ui.input.string('View', {value: this.followedView, placeholder: 'Any view',
      tooltipText: 'Only clicks in the view of this name'});
    const host = ui.box();
    const load = () => {
      this.followedView = (view.value ?? '').trim();
      grok.shell.o = null;
      ui.empty(host);
      host.append(ui.waitBox(async () => {
        try {
          const table = await queries.clicksFollowedByError(filter.date!, filter.groups, this.followedView);
          if (table.rowCount === 0)
            return emptyState('No clicks', 'Choose a longer Date, or clear View');
          table.name = 'Clicks followed by error';
          const grid = DG.Viewer.grid(table, {showRowHeader: false, allowRowSelection: false, allowBlockSelection: false});
          grid.col('element')!.width = 400;
          grid.col('followed_by_error')!.name = 'followed by error ≤ 5 s';
          grid.col('followed_by_error_pct')!.name = 'followed by error, %';
          table.onCurrentRowChanged.subscribe(() => {
            const i = table.currentRowIdx;
            if (i >= 0)
              this.showClickErrors(filter, table.get('element', i), table.get('followed_by_error', i));
          });
          return grid.root;
        }
        catch (e: any) {
          return ui.divText(`Clicks followed by error: ${e?.message ?? e}`, 'd4-viewer-error');
        }
      }));
    };
    view.onChanged.pipe(debounceTime(500)).subscribe(() => load());
    load();
    return ui.divV([ui.div([ui.form([view])], 'ua-toolbar'), host], 'ui-box');
  }

  /** The errors that followed clicks on [element], by signature, each with the timeline of its latest click. */
  showClickErrors(filter: UaFilter, element: string, errors: number): void {
    const acc = DG.Accordion.create();
    acc.addPane('Errors after the click', () => ui.divV([ui.divText(element, 'ua-error-text'), errors === 0 ?
      ui.divText('No error followed these clicks') : ui.wait(async () => {
        const t = await queries.clickErrors(filter.date!, filter.groups, element, this.followedView);
        return rowsTable(t, 'No errors', (r) => [ErrorsView.shortSignature(t.get('signature', r)), t.get('error', r),
          t.get('clicks', r), ui.link('timeline', () => this.openTimeline('action', t.get('action', r)),
            'The latest click')],
        ['signature', 'error', 'clicks', '']);
      })]), true);
    grok.shell.o = acc.root;
  }

  async getClickAnalysisTab(filter: UaFilter): Promise<HTMLDivElement> {
    const table = await queries.getAggregatedClicks(filter.date!);
    table.name = 'Click Analysis';
    this.clicksDf = table;
    const descriptionCol = table.col('description');
    if (!descriptionCol)
      throw new Error('Description column is missing in the Click Analysis table');
    for (let i = 0; i < table.rowCount; i++) {
      const desc = descriptionCol.get(i);
      if (desc.includes('Shell / '))
        descriptionCol.set(i, desc.replace('Shell / ', ''));
    }

    const getPackageName = (path: string) => {
      if (!path)
        return null;
      const parts = path.split(' / ');
      let lastPackage: string | null = null;
      for (const part of parts) {
        const idx = part.indexOf(':');
        if (idx !== -1 && idx < part.length - 1 && part[idx + 1] !== ' ' && isNaN(+part[idx + 1]) &&
          !(part.includes('https:') || part.includes('http:')))
          lastPackage = part.substring(0, idx).trim();
      }
      return lastPackage;
    };

    const sourceCol = table.columns.addNewString('Source');
    const packageNameCol = table.columns.addNewString('Package');
    sourceCol.init((i) => getPackageName(descriptionCol.get(i)) != null ? 'Package' : 'Platform');
    packageNameCol.init((i) => getPackageName(descriptionCol.get(i)));

    const fillLevels = (idx: number, path: string | null) => {
      if (!path)
        return;
      const rawParts = path.split(' / ');
      const finalParts: string[] = [];

      for (const part of rawParts) {
        const colonIndex = part.indexOf(':');
        if (colonIndex !== -1 && colonIndex < part.length - 1 && part[colonIndex + 1] !== ' ')
          finalParts.push(...part.split(':'));
        else if (colonIndex !== -1 && colonIndex < part.length - 1 && part[colonIndex + 1] === ' ') {
          finalParts.push(part.substring(0, colonIndex).trim());
          finalParts.push(part.substring(colonIndex + 1).trim());
        }
        else
          finalParts.push(part);
      }

      for (let i = 0; i < finalParts.length; i++) {
        const colName = `Level ${i + 1}`;
        if (!table.columns.contains(colName))
          table.columns.addNewString(colName);
        table.col(colName)!.set(idx, finalParts[i].trim());
      }
    };

    if (table.rowCount > 0 && descriptionCol) {
      for (let i = 0; i < table.rowCount; i++)
        fillLevels(i, descriptionCol.get(i));
    }

    const grid = DG.Viewer.grid(table, {rowHeight: 20});
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';

    const setGridColWidth = (colName: string, width: number) => {
      const col = grid.col(colName);
      if (col)
        col.width = width;
    };
    setGridColWidth('description', 350);
    setGridColWidth('Level 1', 120);
    setGridColWidth('Level 2', 120);
    setGridColWidth('Level 3', 120);
    grid.sort(['count'], [false]);

    const sourceGridCol = grid.col('Source');
    if (sourceGridCol)
      sourceGridCol.visible = false;
    const packageGridCol = grid.col('Package');
    if (packageGridCol)
      packageGridCol.visible = false;

    const filters = ui.box();
    filters.style.maxWidth = '230px';
    const filtersStyle = {columnNames: ['event_type', 'Source', 'count', 'Package', 'Level 1', 'Level 2', 'Level 3']};
    const filtersRootToAdd = DG.Viewer.filters(table, filtersStyle).root;
    const filtersRootToAddChildren = Array.from(filtersRootToAdd.children);
    for (let i = 0; i < filtersRootToAddChildren.length - 1; i++)
      filtersRootToAddChildren[i].remove();
    filters.append(filtersRootToAdd);
    const treeMapViewer = DG.Viewer.treeMap(table, {
      splitByColumnNames: ['Level 1', 'Level 2', 'Level 3'],
    });
    return ui.splitH([
      filters,
      ui.splitV([
        ui.box(grid.root),
        ui.box(treeMapViewer.root),
      ]),
    ]);
  }
}
