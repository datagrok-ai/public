import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {queries} from '../package-api';
import dayjs from 'dayjs';
import {merge} from 'rxjs';


interface ClickRecord {
  userId: string;
  message: string;
  time: dayjs.Dayjs;
}

export class ClicksView extends UaView {
  expanded: {[key: string]: boolean} = {f: true, l: true};
  tabControl?: DG.TabControl;

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Clicks';
  }

  async initViewers(path?: string): Promise<void> {
    this.root.className = 'grok-view ui-box';
    const tabs: {[key: string]: (() => HTMLElement)} = {
      'Click Analysis': () => this.createWaitElement(this.getClickAnalysisTab, this),
      'Specific clicks': () => this.createWaitElement(this.getSpecificClicksTab, this),
    };
    this.tabControl = ui.tabControl(tabs);
    this.root.appendChild(this.tabControl.root);
  }

  createWaitElement(getElemFunc: (thisObj: ClicksView) => Promise<HTMLElement>, thisObj: ClicksView): HTMLElement {
    const elem: HTMLElement = ui.wait(async () => await getElemFunc(thisObj));
    elem.style.height = '100%';
    return elem;
  }

  async getClickAnalysisTab(thisObj: ClicksView): Promise<HTMLDivElement> {
    const table = await queries.getAggregatedClicks();
    table.name = 'Click Analysis';
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

    const isPackageCol = table.columns.addNewBool('~Is Package');
    const packageNameCol = table.columns.addNewString('~Package Name');
    isPackageCol.init((i) => getPackageName(descriptionCol.get(i)) != null);
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
    setGridColWidth('name', 350);
    setGridColWidth('Level 1', 120);
    setGridColWidth('Level 2', 120);
    setGridColWidth('Level 3', 120);
    grid.sort(['count'], [false]);

    const filters = ui.box();
    filters.style.maxWidth = '230px';
    const filtersStyle = {columnNames: ['~Is Package', '~Package Name', 'Level 1', 'Level 2', 'Level 3', 'count']};
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

  async getSpecificClicksTab(thisObj: ClicksView): Promise<HTMLElement> {
    const PAGE_SIZE = 1000;
    const clickCountBySelector: Map<string, number> = new Map();
    const clickCountByElement: Map<Element, number> = new Map();
    let clickChronology: ClickRecord[] = [];
    let activeZoneOverlays: HTMLElement[] = [];

    const table = DG.DataFrame.create(0, 'Specific Clicks');
    table.columns.addNewString('description');

    const userInput = ui.input.userGroups('User');
    const newUsersChoiceInput = ui.input.choice<DG.User>('User', {items: [], nullable: true});
    newUsersChoiceInput.root.style.display = 'none';
    const newUsersOnlyInput = ui.input.bool('New users only', {value: false});
    const chronologyModeInput = ui.input.bool('Chronology mode', {value: false, tooltipText: 'Show clicks in the order they were made'});
    const highlightZonesInput = ui.input.bool('Highlight zones', {value: false, tooltipText: 'Highlight clicks from different zones'});

    const grid = DG.Viewer.grid(table);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    grid.root.style.position = 'absolute';
    const gridWrapper = ui.div([grid.root], {style: {flex: '1', position: 'relative', overflow: 'hidden'}});
    const lastClickDiv = ui.divV([ui.h3('Last clicked elements:')]);
    lastClickDiv.style.minHeight = '150px';

    const createTable = () => {
      if (table.rowCount > 0)
        table.rows.removeAt(0, table.rowCount);
      const cols = table.columns.toList();
      cols.forEach((c) => table.columns.remove(c.name));
      table.columns.addNewString('description');
      if (chronologyModeInput.value)
        table.columns.addNewDateTime('time');
      else
        table.columns.addNewInt('count');
    };
    createTable();

    const addDataRow = (path: string, options: { count?: number, time?: any } = {}) => {
      const idx = table.rows.addNew().idx;
      table.col('description')!.set(idx, path);
      const countCol = table.columns.byName('count');
      const timeCol = table.columns.byName('time');
      if (options.count != undefined && countCol)
        countCol.set(idx, options.count);
      if (options.time != undefined && timeCol)
        timeCol.set(idx, options.time);
    };

    const increment = (record: ClickRecord) => {
      if (chronologyModeInput.value)
        clickChronology.push(record);
      else {
        const current = clickCountBySelector.get(record.message) || 0;
        clickCountBySelector.set(record.message, current + 1);
      }
    };

    const refresh = () => {
      createTable();
      if (chronologyModeInput.value) {
        const reversedChronology = [...clickChronology].reverse();
        for (const click of reversedChronology)
          addDataRow(click.message, {time: click.time});
        if (table.rowCount > 0)
          grid.sort(['time'], [false]);
      }
      else {
        clickCountBySelector.forEach((count, path) => {
          addDataRow(path, {count});
        });
        if (table.rowCount > 0)
          grid.sort(['count'], [false]);
      }

      if (grid.col('description'))
        grid.col('description')!.width = 500;
      if (grid.col('time')) {
        grid.col('time')!.format = 'yyyy-MM-dd HH:mm:ss';
        grid.col('time')!.width = 150;
      }
      table.fireValuesChanged();
    };

    const addZoneOverlay = (target: Element, color: string) => {
      const rect = target.getBoundingClientRect();
      if (rect.width === 0 || rect.height === 0)
        return;

      const overlay = document.createElement('div');
      overlay.style.position = 'fixed';
      overlay.style.top = `${rect.top}px`;
      overlay.style.left = `${rect.left}px`;
      overlay.style.width = `${rect.width}px`;
      overlay.style.height = `${rect.height}px`;
      overlay.style.zIndex = '2147483640';
      overlay.style.pointerEvents = 'none';
      overlay.style.backgroundColor = color;
      overlay.style.transition = 'opacity 0.2s';

      document.body.appendChild(overlay);
      activeZoneOverlays.push(overlay);
    };

    const clearZoneOverlays = () => {
      activeZoneOverlays.forEach((o) => o.remove());
      activeZoneOverlays = [];
    };

    const updateValues = async () => {
      clickCountBySelector.clear();
      clickCountByElement.clear();
      clickChronology = [];
      highlightZonesInput.value = false;
      highlightZonesInput.enabled = !chronologyModeInput.value;

      const userGroups = userInput.value;
      if (userGroups && userGroups.length > 0) {
        const ids = userGroups.map((ug: any) => `"${ug.user.id}"`).join(',');
        const filter = `eventType.name="click" && session.user.id in (${ids})`;
        const logEventList = await grok.dapi.log
          .filter(filter)
          .list({pageSize: PAGE_SIZE});
        logEventList.forEach((clickEvent: any) => {
          increment({
            userId: clickEvent.session.user.id,
            message: clickEvent.description,
            time: clickEvent.eventTime,
          });
        });
      }
    };

    const setDefaultUser = async () => {
      const currentUser = grok.shell.user;
      const groups = await grok.dapi.groups.getGroupsLookup(currentUser.login);
      const value = groups.find((ug) => ug.personal && ug.user?.login === currentUser.login) || currentUser.group;
      if (value)
        // @ts-ignore
        userInput.value = [value];
    };
    await setDefaultUser();

    DG.debounce(merge(userInput.onChanged, userInput.onInput), 10).subscribe(async () => {
      await updateValues();
      refresh();
    });

    newUsersOnlyInput.onChanged.subscribe(async () => {
      if (newUsersOnlyInput.value) {
        userInput.value = [];
        const allUsers = await grok.dapi.users.order('joined', true).list();
        const cutoff = dayjs().subtract(4, 'day');
        const newUsers = allUsers.filter((u: any) => u.joined && u.joined.isAfter(cutoff));
        newUsersChoiceInput.items = newUsers;
        if (newUsers.length > 0)
          newUsersChoiceInput.value = newUsers[0];
      }
      else {
        newUsersChoiceInput.items = [];
        newUsersChoiceInput.value = null;
        await setDefaultUser();
      }
      userInput.root.style.display = newUsersOnlyInput.value ? 'none' : '';
      newUsersChoiceInput.root.style.display = newUsersOnlyInput.value ? '' : 'none';
    });

    newUsersChoiceInput.onChanged.subscribe(async () => {
      const selectedUser = newUsersChoiceInput.value;
      if (!selectedUser)
        userInput.value = [];
      else {
        const groups = await grok.dapi.groups.getGroupsLookup(selectedUser.login);
        const personalGroup = groups.find((ug) => ug.personal);
        if (personalGroup)
          // @ts-ignore
          userInput.value = [personalGroup];
      }
    });

    chronologyModeInput.onChanged.subscribe(async () => {
      await updateValues();
      createTable();
      refresh();
    });

    highlightZonesInput.onChanged.subscribe(() => {
      if (!highlightZonesInput.value) {
        clearZoneOverlays();
        return;
      }

      const findElements = (segments: string[]): Element[] => {
        let candidates: Element[] = [document.body];
        for (const seg of segments) {
          const nextCandidates: Element[] = [];
          for (const c of candidates) {
            const search = (e: Element) => {
              const desc = DG.ClickUtils.getClickElementDescription(e);
              const dName = e.getAttribute('data-name');
              const dSource = e.getAttribute('data-source');
              const sDesc = desc != null ? DG.ClickUtils.sanitizeCssAttrValue(desc) : null;
              const sName = dName != null ? DG.ClickUtils.sanitizeCssAttrValue(dName) : null;
              const sSource = dSource != null ? DG.ClickUtils.sanitizeCssAttrValue(dSource) : null;

              const check = (s: string | null) => {
                if (!s)
                  return false;
                if (s === seg)
                  return true;
                if (s.endsWith(' ' + seg) || s.endsWith('-' + seg))
                  return true;
                return false;
              };

              if (check(sDesc) || check(sName) || check(sSource))
                nextCandidates.push(e);
              Array.from(e.children).forEach(search);
            };
            search(c);
          }
          candidates = nextCandidates;
          if (candidates.length === 0)
            return [];
        }
        return candidates;
      };

      clickCountByElement.clear();
      clickCountBySelector.forEach((count, description) => {
        const cleanDesc = description.replace(' (right button)', '').replace(' (middle button)', '');
        const segments = cleanDesc.split(' / ');
        if (segments.length === 1 && segments[0] === 'Shell')
          return;
        const elements = findElements(segments);
        if (!elements || elements.length === 0)
          return;
        elements.forEach((el) => {
          const current = clickCountByElement.get(el) || 0;
          clickCountByElement.set(el, current + count);
        });
      });

      const counts = Array.from(clickCountByElement.values());
      let minCount = 1, maxCount = 1;
      if (counts.length > 0) {
        minCount = Math.min(...counts);
        maxCount = Math.max(...counts);
      }

      clearZoneOverlays();
      clickCountByElement.forEach((totalClicks, element) => {
        const colorScheme = [DG.Color.lightLightGray, DG.Color.barChart, DG.Color.success, DG.Color.darkGreen];
        const finalClicks = Math.min(Math.max(totalClicks, minCount), maxCount);
        const below = (colorScheme[0] & 0x00FFFFFF) | (75 << 24);
        const above = (colorScheme[colorScheme.length - 1] & 0x00FFFFFF) | (75 << 24);
        const colorToUse = DG.Color.scaleColor(finalClicks, minCount, maxCount, 75, colorScheme, below, above);
        const selectorColor = `rgba(${DG.Color.r(colorToUse)}, ${DG.Color.g(colorToUse)}, ${DG.Color.b(colorToUse)}, ${DG.Color.a(colorToUse) / 255})`;
        addZoneOverlay(element, selectorColor);
      });
    });

    thisObj.sub(grok.events.onLog.subscribe((msg: DG.LogMessage) => {
      if (msg.level !== 'usage' || msg.auditType !== 'click')
        return;
      const userGroups = userInput.value;
      const currentUserId = grok.shell.user.id;
      if (userGroups && userGroups.length > 0 && userGroups.some((ug: DG.User) => ug./*user.*/id === currentUserId)) {
        const record: ClickRecord = {
          userId: currentUserId,
          message: msg.message,
          time: DG.toJs(msg.time),
        };
        increment(record);

        if (chronologyModeInput.value) {
          addDataRow(record.message, {time: record.time});
          grid.sort(['time'], [false]);
        }
        else {
          let idx = -1;
          const nameCol = table.col('name');
          const countCol = table.col('count');
          if (nameCol) {
            for(let i = 0; i < table.rowCount; i++) {
              if (nameCol.get(i) === record.message) {
                idx = i;
                break;
              }
            }
          }
          if (idx === -1)
            addDataRow(record.message, {count: 1});
          else if (countCol)
            countCol.set(idx, (countCol.get(idx) || 0) + 1);
          grid.sort(['count'], [false]);
        }
        table.fireValuesChanged();
      }

      if (lastClickDiv.children.length > 5)
        lastClickDiv.lastChild?.remove();
      const entry = ui.divText(`${msg.message}`);
      if (lastClickDiv.children.length > 0)
        lastClickDiv.insertBefore(entry, lastClickDiv.children[1]);
      else
        lastClickDiv.appendChild(entry);
    }));

    return ui.divV([
      userInput.root,
      newUsersChoiceInput.root,
      newUsersOnlyInput.root,
      chronologyModeInput.root,
      highlightZonesInput.root,
      gridWrapper,
      lastClickDiv,
    ], {classes: 'grok-inspector', style: {marginLeft: '12px'}});
  }
}
