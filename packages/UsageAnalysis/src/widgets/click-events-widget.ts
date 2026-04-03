import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {merge} from 'rxjs';
import dayjs from 'dayjs';


interface ClickRecord {
  userId: string;
  message: string;
  time: dayjs.Dayjs;
}

export class ClickEventsWidget extends DG.Widget {
  readonly ready: Promise<void>;
  pickModeCleanup: (() => void) | null = null;

  constructor() {
    super(ui.div([], {style: {height: '100%'}}));
    this.ready = this._buildUI().then((content) => {
      this.root.replaceChildren(content);
    });
  }

  private async _buildUI(): Promise<HTMLElement> {
    const PAGE_SIZE = 1000;
    const clickCountBySelector: Map<string, number> = new Map();
    const clickCountByElement: Map<Element, number> = new Map();
    let clickChronology: ClickRecord[] = [];
    let activeZoneOverlays: HTMLElement[] = [];

    const table = DG.DataFrame.create(0, 'Click Events');
    table.columns.addNewString('description');

    const userInput = ui.input.userGroups('User');
    const newUsersChoiceInput = ui.input.choice<DG.User>('User', {items: [], nullable: true});
    newUsersChoiceInput.root.style.display = 'none';
    const newUsersOnlyInput = ui.input.bool('New users only', {value: false});
    const chronologyModeInput = ui.input.bool('Latest clicks', {value: false, tooltipText: 'Show individual click events sorted by time, newest first'});
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
        ui.setUpdateIndicator(gridWrapper, true);
        try {
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
        } finally {
          ui.setUpdateIndicator(gridWrapper, false);
        }
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

    this.subs.push(DG.debounce(merge(userInput.onChanged, userInput.onInput), 10).subscribe(async () => {
      await updateValues();
      refresh();
    }));

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

    this.subs.push(grok.events.onLog.subscribe((msg: DG.LogMessage) => {
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
            for (let i = 0; i < table.rowCount; i++) {
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

    const pickerIcon = ui.iconFA('crosshairs', (e) => {
      e.preventDefault();
      this.togglePicker(e.currentTarget as HTMLElement);
    }, 'Pick element to analyze usage');
    pickerIcon.style.cursor = 'pointer';

    const form = ui.form([userInput, newUsersChoiceInput, newUsersOnlyInput, chronologyModeInput, highlightZonesInput]);
    form.style.flex = '0 0 auto';

    return ui.divV([
      ui.divH([pickerIcon], {style: {padding: '2px 0 2px 4px'}}),
      form,
      gridWrapper,
      lastClickDiv,
    ], {classes: 'grok-inspector', style: {marginLeft: '12px'}});
  }

  async showElementStats(selector: string) {
    document.body.style.cursor = 'wait';
    try {
      const filter = `eventType.name="click" && description like "%${selector}%"`;
      const logs = await grok.dapi.log.filter(filter).list();
      const userIds: Set<string> = new Set<string>();
      for (const log of logs)
        if (log.session && log.session instanceof DG.UserSession && log.session.user && log.session.user.id)
          userIds.add(log.session.user.id);
      const idToName: Map<string, string> = new Map<string, string>();
      await Promise.all(Array.from(userIds).map(async (id) => {
        try {
          if (!id)
            return;
          const user = await grok.dapi.users.find(id);
          idToName.set(id, user.friendlyName ?? user.firstName ?? user.login ?? 'Unknown');
        } catch (e) {
          idToName.set(id!, `Unknown (${id})`);
        }
      }));

      const userCounts: {[key: string]: number} = {};
      for (const log of logs) {
        if (log.session && log.session instanceof DG.UserSession) {
          const id = log.session.user?.id;
          const displayName = (id && idToName.has(id)) ? idToName.get(id)! : 'Unknown';
          userCounts[displayName] = (userCounts[displayName] || 0) + 1;
        }
      }
      const statsTable = DG.DataFrame.fromColumns([
        DG.Column.fromStrings('User', Object.keys(userCounts)),
        DG.Column.fromList(DG.TYPE.INT, 'Clicks', Object.values(userCounts))]);

      ui.dialog('Element Usage')
        .add(ui.h1(selector))
        .add(DG.Viewer.grid(statsTable).root)
        .show();
    } finally {
      document.body.style.cursor = 'default';
    }
  }

  togglePicker(iconRoot: HTMLElement) {
    let highlighter: HTMLDivElement | null = null;
    const removeHighlight = () => {
      highlighter?.remove();
      highlighter = null;
    };

    const drawHighlight = (target: HTMLElement) => {
      removeHighlight();
      const rect = target.getBoundingClientRect();
      highlighter = document.createElement('div');
      Object.assign(highlighter.style, {
        position: 'fixed',
        top: `${rect.top}px`,
        left: `${rect.left}px`,
        width: `${rect.width}px`,
        height: `${rect.height}px`,
        zIndex: '2147483647',
        pointerEvents: 'none',
        boxSizing: 'border-box',
        border: '2px solid rgba(50, 150, 255, 1)',
        backgroundColor: 'rgba(50, 150, 255, 0.3)'
      });
      document.body.appendChild(highlighter);
    };

    const disable = () => {
      if (this.pickModeCleanup) {
        this.pickModeCleanup();
        this.pickModeCleanup = null;
      }
      document.body.style.cursor = 'default';
      iconRoot.style.color = '';
      removeHighlight();
    };

    if (this.pickModeCleanup)
      disable();
    else {
      document.body.style.cursor = 'crosshair';
      iconRoot.style.color = 'var(--blue-1)';

      const handleMouseDown = (e: MouseEvent) => {
        const target = e.target as HTMLElement;
        if (target === iconRoot || iconRoot.contains(target))
          return;
        e.preventDefault();
        e.stopPropagation();
        disable();

        const handleClick = (eClick: MouseEvent) => {
          eClick.preventDefault();
          eClick.stopPropagation();
          document.removeEventListener('click', handleClick, true);
        };
        document.addEventListener('click', handleClick, true);

        setTimeout(() => {
          document.removeEventListener('click', handleClick, true);
        }, 300);

        const menu = DG.Menu.popup();
        menu.onClose.subscribe(() => removeHighlight());
        let currentTarget: HTMLElement | null = target;
        const itemTitles = new Set<string>();
        while (currentTarget && currentTarget !== document.body) {
          const logName = DG.ClickUtils.getElementLoggingName(currentTarget!);
          if (logName && logName !== '' && logName !== 'Shell') {
            const fullPath = DG.ClickUtils.getFullPath(currentTarget!);
            const tagName = currentTarget!.tagName.toLowerCase();
            const elem = currentTarget;
            const itemTitle = `${logName} (${tagName})`;
            if (!itemTitles.has(itemTitle)) {
              itemTitles.add(itemTitle);
              menu.item(itemTitle, () => this.showElementStats(fullPath), null, {onMouseEnter: () => drawHighlight(elem!),
                onMouseLeave: () => removeHighlight()});
            }
          }
          currentTarget = currentTarget!.parentElement;
        }
        menu.show({causedBy: e});
      };
      document.addEventListener('mousedown', handleMouseDown, true);
      this.pickModeCleanup = () => {
        document.removeEventListener('mousedown', handleMouseDown, true);
        document.body.style.cursor = 'default';
        iconRoot.style.color = '';
      };
    }
  }
}
