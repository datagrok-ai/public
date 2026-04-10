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

const MAX_LOADED_EVENTS = 10_000;

export class ClickEventsWidget extends DG.Widget {
  readonly ready: Promise<void>;
  pickModeCleanup: (() => void) | null = null;

  constructor() {
    super(ui.div([], {style: {height: '100%', width: '100%'}}));
    this.ready = this._buildUI().then((content) => {
      this.root.replaceChildren(content);
    });
  }

  private async _buildUI(): Promise<HTMLElement> {
    const clickCountBySelector: Map<string, number> = new Map();
    const clickCountByElement: Map<Element, number> = new Map();
    let clickChronology: ClickRecord[] = [];
    let activeZoneOverlays: HTMLElement[] = [];

    const table = DG.DataFrame.create(0, 'Click Events');
    table.columns.addNewString('description');

    const userInput = ui.input.user('User', {tooltipText: 'User whose click activity to analyze'});
    let newUsersCache: DG.User[] = [];
    const newUsersTagsInput = ui.input.tags<DG.User>('User', {
      tooltipText: 'Select specific new users to analyze; leave empty to show no data',
      getSuggestions: async (text: string) => {
        const lower = text.toLowerCase();
        return newUsersCache.filter((u) => {
          const name = (u.friendlyName ?? u.login ?? '').toLowerCase();
          return !lower || name.includes(lower);
        });
      },
    });
    newUsersTagsInput.root.style.display = 'none';
    const newUsersOnlyInput = ui.input.bool('New users only', {value: false, tooltipText: 'Filter to users who joined in the last 4 days'});
    const eventTypeInput = ui.input.choice('Event type', {
      items: ['Clicks', 'Menu clicks', 'Dialogs', 'Inputs', 'Commands', 'All'],
      value: 'Clicks',
      tooltipText: 'Type of platform interactions to analyze',
    });
    const dateRangeInput = ui.input.choice('Date range', {
      items: ['Today', 'Last 7 days', 'Last 30 days', 'All time'],
      value: 'Last 30 days',
      tooltipText: 'Time window for loading historical data',
    });
    const chronologyModeInput = ui.input.bool('Latest events', {value: false, tooltipText: 'Show individual events sorted by time, newest first'});
    const highlightZonesInput = ui.input.bool('Highlight zones', {value: false, tooltipText: 'Overlay color-coded highlights on UI elements showing their click frequency'});

    const eventTypeFilters: Record<string, string> = {
      'Clicks': 'eventType.name="click"',
      'Menu clicks': 'eventType.name="menu click"',
      'Dialogs': 'eventType.name in ("dialog show","dialog ok","dialog close")',
      'Inputs': 'eventType.name="input"',
      'Commands': 'eventType.name="command"',
      'All': 'eventType.name in ("click","menu click","dialog show","dialog ok","dialog close","input","command","navigate")',
    };
    const eventTypeAuditTypes: Record<string, string[]> = {
      'Clicks': ['click'],
      'Menu clicks': ['menu click'],
      'Dialogs': ['dialog show', 'dialog ok', 'dialog close'],
      'Inputs': ['input'],
      'Commands': ['command'],
      'All': ['click', 'menu click', 'dialog show', 'dialog ok', 'dialog close', 'input', 'command', 'navigate'],
    };

    const grid = DG.Viewer.grid(table);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    grid.root.style.position = 'absolute';
    const gridWrapper = ui.div([grid.root], {classes: 'clicks-grid-wrapper'});

    const eventCountLabel = ui.div([], {classes: 'clicks-event-count'});
    eventCountLabel.style.display = 'none';

    const activityStreamLabel = ui.div(['Activity stream'], {classes: 'clicks-activity-label'});
    const activityStreamEmpty = ui.div(['Interact with the platform to see events here'], {classes: 'clicks-activity-empty'});
    const activityStreamDiv = ui.divV([activityStreamLabel, activityStreamEmpty]);
    activityStreamDiv.classList.add('clicks-activity-stream');

    const createTable = () => {
      if (table.rowCount > 0)
        table.rows.removeAt(0, table.rowCount);
      const cols = table.columns.toList();
      cols.forEach((c) => table.columns.remove(c.name));
      table.columns.addNewString('element');
      if (chronologyModeInput.value)
        table.columns.addNewDateTime('time');
      else
        table.columns.addNewInt('count');
    };
    createTable();

    const addDataRow = (path: string, options: { count?: number, time?: dayjs.Dayjs } = {}) => {
      const idx = table.rows.addNew().idx;
      table.col('element')!.set(idx, path);
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
        clickCountBySelector.forEach((count, path) => addDataRow(path, {count}));
        if (table.rowCount > 0)
          grid.sort(['count'], [false]);
      }

      grid.columns.rowHeader!.width = grid.dataFrame.rowCount.toString().length * 8 + 12;
      if (grid.col('element'))
        grid.col('element')!.width = 260;
      if (grid.col('time')) {
        grid.col('time')!.format = 'yyyy-MM-dd HH:mm:ss';
        grid.col('time')!.width = 140;
      }
      table.fireValuesChanged();

      const count = chronologyModeInput.value ? clickChronology.length : clickCountBySelector.size;
      eventCountLabel.textContent = `${count} event${count !== 1 ? 's' : ''}`;
      eventCountLabel.style.display = count > 0 ? '' : 'none';
    };

    const addZoneOverlay = (target: Element, color: string) => {
      const rect = target.getBoundingClientRect();
      if (rect.width === 0 || rect.height === 0)
        return;

      const overlay = ui.div([], {classes: 'clicks-zone-overlay'});
      overlay.style.top = `${rect.top}px`;
      overlay.style.left = `${rect.left}px`;
      overlay.style.width = `${rect.width}px`;
      overlay.style.height = `${rect.height}px`;
      overlay.style.backgroundColor = color;

      document.body.appendChild(overlay);
      activeZoneOverlays.push(overlay);
    };

    const clearZoneOverlays = () => {
      activeZoneOverlays.forEach((o) => o.remove());
      activeZoneOverlays = [];
    };

    const updateHighlightZonesState = () => {
      const available = !chronologyModeInput.value && eventTypeInput.value === 'Clicks';
      highlightZonesInput.enabled = available;
      highlightZonesInput.setTooltip(available ? 'Highlight elements by click frequency' : chronologyModeInput.value ?
        'Not available in Latest events mode' : 'Only available for the Clicks event type');
    };

    const buildDateFilter = (): string => {
      const fmt = (d: dayjs.Dayjs) => d.format('YYYY-MM-DDTHH:mm:ss');
      const now = dayjs();
      switch (dateRangeInput.value) {
        case 'Today': return `eventTime >= "${fmt(now.startOf('day'))}"`;
        case 'Last 7 days': return `eventTime >= "${fmt(now.subtract(7, 'day').startOf('day'))}"`;
        case 'Last 30 days': return `eventTime >= "${fmt(now.subtract(30, 'day').startOf('day'))}"`;
        default: return '';
      }
    };

    const getActiveUsers = (): Array<DG.User | DG.Group> => {
      if (newUsersOnlyInput.value)
        return (newUsersTagsInput.value as unknown as DG.User[] | null) ?? [];
      return (userInput.value as unknown as DG.Group[] | null) ?? [];
    };

    const toUserId = (u: DG.User | DG.Group): string => {
      const groupUser = (u as DG.Group).user as unknown as DG.User | null;
      return groupUser != null ? groupUser.id : u.id;
    };

    const updateValues = async () => {
      clickCountBySelector.clear();
      clickCountByElement.clear();
      clickChronology = [];
      highlightZonesInput.value = false;
      updateHighlightZonesState();

      const activeUsers = getActiveUsers();
      const ids = activeUsers.map((u) => toUserId(u)).filter((id) => id.length > 0).map((id) => `"${id}"`).join(',');
      if (ids) {
        ui.setUpdateIndicator(gridWrapper, true);
        try {
          const eventFilter = eventTypeFilters[eventTypeInput.value ?? 'Clicks'];
          const dateFilter = buildDateFilter();
          const filter = `${eventFilter} && session.user.id in (${ids})${dateFilter ? ' && ' + dateFilter : ''}`;
          const logEventList = await grok.dapi.log
            .filter(filter)
            .order('eventTime', false)
            .by(MAX_LOADED_EVENTS)
            .list();
          let minTime: dayjs.Dayjs | null = null;
          let maxTime: dayjs.Dayjs | null = null;
          for (const clickEvent of logEventList) {
            if (!clickEvent.description?.trim())
              continue;
            if (!(clickEvent.session instanceof DG.UserSession))
              continue;
            increment({
              userId: clickEvent.session.user.id,
              message: clickEvent.description,
              time: clickEvent.eventTime,
            });
            if (minTime === null || clickEvent.eventTime.isBefore(minTime))
              minTime = clickEvent.eventTime;
            if (maxTime === null || clickEvent.eventTime.isAfter(maxTime))
              maxTime = clickEvent.eventTime;
          }
          if (minTime !== null && maxTime !== null) {
            const fmt = (d: dayjs.Dayjs) => d.format('MMM D, YYYY');
            const rangeStr = minTime.isSame(maxTime, 'day') ? fmt(minTime) : `${fmt(minTime)} – ${fmt(maxTime)}`;
            ui.tooltip.bind(eventCountLabel, rangeStr);
          }
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
        const allUsers = await grok.dapi.users.order('joined', true).list();
        const cutoff = dayjs().subtract(4, 'day');
        newUsersCache = allUsers.filter((u) => u.joined && u.joined.isAfter(cutoff));
      }
      else {
        newUsersCache = [];
        await setDefaultUser();
      }
      userInput.root.style.display = newUsersOnlyInput.value ? 'none' : '';
      newUsersTagsInput.root.style.display = newUsersOnlyInput.value ? '' : 'none';
      await updateValues();
      refresh();
    });

    this.subs.push(DG.debounce(merge(newUsersTagsInput.onChanged, newUsersTagsInput.onInput), 10).subscribe(async () => {
      await updateValues();
      refresh();
    }));

    chronologyModeInput.onChanged.subscribe(async () => {
      updateHighlightZonesState();
      activityStreamDiv.style.display = chronologyModeInput.value ? 'none' : '';
      await updateValues();
      createTable();
      refresh();
    });

    eventTypeInput.onChanged.subscribe(async () => {
      if (highlightZonesInput.value && eventTypeInput.value !== 'Clicks') {
        highlightZonesInput.value = false;
        clearZoneOverlays();
      }
      updateHighlightZonesState();
      await updateValues();
      refresh();
    });

    dateRangeInput.onChanged.subscribe(async () => {
      await updateValues();
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
      if (msg.level !== 'usage')
        return;
      if (!msg.message || !msg.message.trim())
        return;
      const auditTypes = eventTypeAuditTypes[eventTypeInput.value ?? 'Clicks'];
      if (msg.auditType && !auditTypes.includes(msg.auditType))
        return;
      const currentUserId = grok.shell.user.id;
      if (getActiveUsers().some((u) => toUserId(u) === currentUserId)) {
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
          const elemCol = table.col('element');
          const countCol = table.col('count');
          if (elemCol) {
            for (let i = 0; i < table.rowCount; i++) {
              if (elemCol.get(i) === record.message) {
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
        const count = chronologyModeInput.value ? clickChronology.length : clickCountBySelector.size;
        eventCountLabel.textContent = `${count} event${count !== 1 ? 's' : ''}`;
        eventCountLabel.style.display = '';
      }

      const entry = ui.div([msg.message], {classes: 'clicks-activity-entry'});
      activityStreamEmpty.style.display = 'none';
      if (activityStreamDiv.children.length > 1)
        activityStreamDiv.insertBefore(entry, activityStreamDiv.children[1]);
      else
        activityStreamDiv.appendChild(entry);
      while (activityStreamDiv.children.length > 6)
        activityStreamDiv.lastChild?.remove();
    }));

    const pickerIcon = ui.iconFA('crosshairs', undefined, 'Pick element to analyze usage');
    const pickerRow = ui.divH([pickerIcon, ui.span(['Pick element'])], {classes: 'clicks-picker-row'});
    pickerRow.addEventListener('click', (e) => {
      e.preventDefault();
      this.togglePicker(pickerIcon);
    });

    const sectionLabel = (text: string) => ui.div([text], {classes: 'clicks-section-label'});

    const filtersForm = ui.form([userInput, newUsersTagsInput, newUsersOnlyInput, eventTypeInput, dateRangeInput]);
    filtersForm.classList.add('clicks-filters-form');

    const viewRow = ui.divH([chronologyModeInput.root, highlightZonesInput.root], {classes: 'clicks-view-row'});

    updateHighlightZonesState();

    return ui.divV([
      sectionLabel('Filters'),
      filtersForm,
      sectionLabel('View'),
      viewRow,
      pickerRow,
      eventCountLabel,
      gridWrapper,
      activityStreamDiv,
    ], {classes: 'grok-inspector clicks-widget'});
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
      highlighter = ui.div([], {classes: 'clicks-pick-highlight'});
      Object.assign(highlighter.style, {
        top: `${rect.top}px`,
        left: `${rect.left}px`,
        width: `${rect.width}px`,
        height: `${rect.height}px`,
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
