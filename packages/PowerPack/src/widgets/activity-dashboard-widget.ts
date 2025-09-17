import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';
import {queries} from '../package-api';


enum SpotlightTabNames {
  ACTION_REQUIRED = 'Action required',
  SHARED_WITH_ME = 'Shared with me',
  RECENT = 'Recent',
  ADMIN = 'Admin',
}

export class ActivityDashboardWidget extends DG.Widget {
  static RECENT_TIME_DAYS = 2;
  static SPOTLIGHT_ITEMS_LENGTH: 10;

  static sharedEntityRegex: RegExp = /<span>#{x\.([a-f0-9-]+)\.".*?"}<\/span>/g;
  keywordsToIgnore: string[] = ['"data"', '"data w/', '"fitted data"', '"reduced data"', '"reduceddata"', '"sizingparams"', '"primaryfilter', '"trimmed data"', '"input data from imported file"', '"dataanalysisdf', '"clean data"'/*, 'Ran <span>#{x.'*/];
  tipsOfTheDay: string[] = ['Double-click a column header to sort the data.', 'Drag and drop a CSV or Excel file into Datagrok to load it as a table.', 'Hover over a column name to see its statistics in the tooltip.', 'Hold Shift and click multiple column headers to select several columns.', 'Use Ctrl+F in a table to search across all values.', 'Resize columns by dragging the separator in the header.', 'Use Layouts to save and restore your favorite viewer arrangements.', 'Use “Add New Column” to create calculated columns with formulas.', 'Right-click a column to see context actions available for it.', 'Try the “Correlation Plot” to quickly see relationships between numeric columns.'];
  cutoffDate: dayjs.Dayjs = dayjs().subtract(ActivityDashboardWidget.RECENT_TIME_DAYS, 'day');

  favoritesEvents: DG.LogEvent[] = [];

  recentNotifications: DG.UserNotification[] = [];
  recentUserActivity: DG.LogEvent[] = [];

  sharedNotifications: DG.UserNotification[] = [];
  sharedUsers: DG.User[] = [];
  sharedWithMe: DG.Entity[] = [];
  recentEntities: DG.Entity[] = [];
  recentEntityTimes: dayjs.Dayjs[] = [];

  tabControl?: DG.TabControl;
  spotlightRoot?: HTMLElement;
  favoritesListRoot?: HTMLElement;
  notificationsListRoot?: HTMLElement;
  activityListRoot?: HTMLElement;

  actionRequiredRoot?: HTMLDivElement | null = null;
  sharedWithMeRoot?: HTMLDivElement | null = null;
  spacesRoot?: HTMLDivElement | null = null;
  recentItemsRoot?: HTMLDivElement | null = null;
  adminRoot?: HTMLDivElement | null = null;
  subwidgetsAdded: Map<string, HTMLDivElement | null> = new Map<string, HTMLDivElement>();
  subwidgetsAmount: number = 0;
  reportAmount?: number = 0;

  constructor() {
    super(ui.panel([], 'power-pack-activity-widget-content'));
    this.buildTabbedUI();
  }


  async buildTabbedUI(): Promise<void> {
    const tabs: {[key: string]: (() => HTMLElement)} = {
      'Spotlight': () => ui.wait(async () => await this.getSpotlightTab()),
      'Favorites': () => ui.wait(async () => await this.getFavoritesTab()),
      'Notifications': () => ui.wait(async () => await this.getNotificationsTab()),
      'My Activity': () => ui.wait(async () => await this.getActivityTab()),
    };

    this.tabControl = ui.tabControl(tabs, true);
    this.tabControl.onTabChanged.subscribe((_) => this.cleanLists());
    this.tabControl.root.style.height = '100%';
    this.tabControl.root.style.width = '100%';
    (this.tabControl.root.querySelector('.ui-box.d4-tab-content') as HTMLElement)!.style.overflow = 'scroll';

    this.root.appendChild(this.tabControl.root);
  }

  async initSpotlightData(): Promise<void> {
    console.time('ActivityDashboardWidget.initSpotlightData');
    const notificationsDataSource: DG.NotificationsDataSource = grok.dapi.users.notifications.forCurrentUser()
      .by(ActivityDashboardWidget.SPOTLIGHT_ITEMS_LENGTH) as DG.NotificationsDataSource;
    console.time('ActivityDashboardWidget.notificationsAndMostRecentEntities');
    const [notifications, mostRecentEntitiesDf]: [DG.UserNotification[], DG.DataFrame] = await Promise.all([notificationsDataSource.list({pageSize: 20}),
      queries.mostRecentEntities(DG.User.current().id)]);
    console.timeEnd('ActivityDashboardWidget.notificationsAndMostRecentEntities');

    this.recentNotifications = notifications
      .filter((n) => !n.isRead || n.createdAt?.isAfter(this.cutoffDate))
      .sort((a, b) => b.createdAt?.diff(a.createdAt) ?? 0);
    const sharedCandidates = this.recentNotifications.filter((n) => n.text.startsWith('<span>#') &&
      n.text.includes('</span> shared') && n.text.includes(': <span>')).slice(0, ActivityDashboardWidget.SPOTLIGHT_ITEMS_LENGTH);
    const sharedUserIds: string[] = [];
    const sharedEntityIds: string[] = [];
    sharedCandidates.forEach((n) => {
      ActivityDashboardWidget.sharedEntityRegex.lastIndex = 0;
      const matches = [...n.text.matchAll(ActivityDashboardWidget.sharedEntityRegex)];
      if (matches.length > 1) {
        sharedUserIds.push(matches[0][1]);
        sharedEntityIds.push(matches[1][1]);
      }
    });
    this.sharedNotifications = sharedCandidates;

    const recentEntityIdCol = mostRecentEntitiesDf.col('id')!;
    const lastEventTimeCol = mostRecentEntitiesDf.col('last_event_time')!;
    const recentEntityIds = new Array(recentEntityIdCol.length);
    for (let i = 0; i < recentEntityIdCol.length; i++)
      recentEntityIds[i] = recentEntityIdCol.get(i);

    const uniqueIds = Array.from(new Set<string>([...sharedUserIds, ...sharedEntityIds, ...recentEntityIds]));
    console.time('ActivityDashboardWidget.getEntitiesByIds');
    const allEntities = await grok.dapi.getEntities(uniqueIds);
    console.timeEnd('ActivityDashboardWidget.getEntitiesByIds');

    const byIdMap = new Map(allEntities.map((e) => [e.id, e]));

    this.sharedUsers = [];
    this.sharedWithMe = [];
    for (let i = 0; i < sharedUserIds.length; i++) {
      const user = byIdMap.get(sharedUserIds[i]);
      if (user instanceof DG.User)
        this.sharedUsers.push(user);
    }
    for (let i = 0; i < sharedEntityIds.length; i++) {
      const ent = byIdMap.get(sharedEntityIds[i]);
      if (!!ent && !(ent instanceof DG.User))
        this.sharedWithMe.push(ent);
    }
    this.recentEntities.length = 0;
    this.recentEntityTimes.length = 0;
    for (let i = 0; i < recentEntityIds.length; i++) {
      const ent = byIdMap.get(recentEntityIds[i]);
      if (!ent)
        continue;
      if (!(ent instanceof DG.FuncCall || ent instanceof DG.Group || ent instanceof DG.User || ent instanceof DG.Package ||
        ent instanceof DG.UserReport || ent instanceof DG.TableInfo || (ent instanceof DG.Func && !(ent instanceof DG.Script ||
        ent instanceof DG.DataQuery || ent instanceof DG.DataJob)) || ent instanceof DG.ViewInfo || ent instanceof DG.DataConnection || ent == null ||
        //@ts-ignore
        (ent instanceof DG.Project && (!ent.isDashboard || ent.isPackage)) || (ent.hasOwnProperty('npmScope') && ent['npmScope'] == 'datagrok'))) {
        this.recentEntities.push(ent);
        this.recentEntityTimes.push(lastEventTimeCol.get(i));
      }
    }
    console.timeEnd('ActivityDashboardWidget.initSpotlightData');
  }

  async getSpotlightTab(): Promise<HTMLElement> {
    console.time('ActivityDashboardWidget.buildSpotlightTab');
    await this.initSpotlightData();

    const createSection = (title: SpotlightTabNames, items: DG.Entity[] | DG.UserNotification[], icon: HTMLElement) => {
      let usedItems: DG.Entity[] | DG.UserNotification[] = [];
      if (title === SpotlightTabNames.ACTION_REQUIRED) {
        this.reportAmount = (items as DG.UserNotification[])
          .filter((item) => item.text.toLowerCase().includes('you were assigned')).length;
        let reportPresent = false;
        if (this.reportAmount > 0) {
          for (let i = 0; i < items.length; i++) {
            const item = items[i] as DG.UserNotification;
            if (item.text.toLowerCase().includes('you were assigned')) {
              if (!reportPresent) {
                usedItems.push(item as (DG.UserNotification & DG.Entity));
                reportPresent = true;
              }
            } else
              usedItems.push(item as (DG.UserNotification & DG.Entity));
          }
        }
      } else
        usedItems = items;
      const itemsToShow = usedItems.slice(0, ActivityDashboardWidget.SPOTLIGHT_ITEMS_LENGTH);
      const list = ui.list(itemsToShow);
      list.classList.add('power-pack-activity-widget-subwidget-list-content');
      for (let i = 0; i < itemsToShow.length; i++) {
        const item = itemsToShow[i];
        const listChild = list.children[i];
        const sharedNotification = this.sharedNotifications[i];
        if (item instanceof DG.Project || item instanceof DG.Notebook || item instanceof DG.Model || item instanceof DG.ViewLayout) {
          const parentElem = listChild.firstChild;
          const iconToPaste = item instanceof DG.Project ? ui.iconSvg('project') : item instanceof DG.Notebook ? ui.iconImage('notebook', '/images/entities/jupyter.png') :
            item instanceof DG.Model ? ui.iconSvg('model') : ui.iconSvg('view-layout');
          if (parentElem && parentElem.firstChild)
            parentElem.insertBefore(iconToPaste, parentElem.firstChild);
        } else if (item instanceof DG.UserNotification) {
          if (item.createdAt) {
            const timestamp = ui.time.shortTimestamp(item.createdAt);
            timestamp.style.top = '5px';
            listChild.prepend(timestamp);
            const timeChild = listChild.querySelector('.d4-time');
            if (timeChild)
              timeChild.remove();
          }
        } else if (item instanceof DG.LogEvent && item.description.toLowerCase().includes('published version'))
          listChild.querySelector('.d4-markup')!.classList.add('power-pack-activity-widget-spotlight-column-admin-row');

        if (title === SpotlightTabNames.SHARED_WITH_ME && sharedNotification.text) {
          const icon = ui.iconFA('clock');
          ui.tooltip.bind(icon, () => {
            const user = this.sharedUsers[i];
            const entity = this.sharedWithMe[i];
            const inputText = ui.span([]);
            inputText.innerHTML = sharedNotification.text;
            const userPicture = user.picture as HTMLElement;
            const entityPicture = item instanceof DG.Project ? ui.iconSvg('project') : item instanceof DG.Notebook ? ui.iconImage('notebook', '/images/entities/jupyter.png') :
              item instanceof DG.Model ? ui.iconSvg('model') : ui.iconSvg('view-layout');
            //@ts-ignore
            return this.replaceMarkupSpans(inputText, [userPicture, entityPicture], [user.friendlyName, entity.friendlyName], sharedNotification.createdAt);
          });
          icon.style.top = '5px';
          listChild.prepend(icon);
        } else if (title === SpotlightTabNames.RECENT) {
          const timestamp = ui.time.shortTimestamp(this.recentEntityTimes[i]);
          timestamp.style.top = '5px';
          listChild.prepend(timestamp);
        }
      }
      this.subwidgetsAmount++;
      return ui.divV([ui.h3(ui.span([icon, ui.span([` ${title}`])]), 'power-pack-activity-widget-spotlight-column-header'),
        list], title === SpotlightTabNames.ACTION_REQUIRED ? 'power-pack-activity-widget-spotlight-column-action-required' : 'power-pack-activity-widget-spotlight-column');
    };

    let root = ui.divH([], 'power-pack-activity-widget-spotlight-root');

    // API changes needed
    // if (!(DG.User.current().joined > dayjs().subtract(7, 'day')) && false) {
    const actionRequired = this.recentNotifications.filter((n) => {
      const text = n.text.toLowerCase();
      return text.includes('you were assigned') || text.includes('requested a membership');
    });
    this.actionRequiredRoot = actionRequired.length > 0 ? createSection(SpotlightTabNames.ACTION_REQUIRED, actionRequired, ui.iconFA('exclamation-circle')) : null;
    if (this.actionRequiredRoot)
      root.appendChild(this.actionRequiredRoot);
    this.subwidgetsAdded.set(SpotlightTabNames.ACTION_REQUIRED, this.actionRequiredRoot);

    this.sharedWithMeRoot = this.sharedWithMe.length > 0 ? createSection(SpotlightTabNames.SHARED_WITH_ME, this.sharedWithMe, ui.iconFA('inbox')) : null;
    if (this.sharedWithMeRoot)
      root.appendChild(this.sharedWithMeRoot);
    this.subwidgetsAdded.set(SpotlightTabNames.SHARED_WITH_ME, this.sharedWithMeRoot);

    this.spacesRoot = null;
    if (this.spacesRoot)
      root.appendChild(this.spacesRoot);
    // this.subwidgetsAdded.set(SpotlightTabNames.SPACES, this.spacesRoot);

    this.recentItemsRoot = this.recentEntities.length > 0 ? createSection(SpotlightTabNames.RECENT, this.recentEntities, ui.iconFA('history')) : null;
    if (this.recentItemsRoot)
      root.appendChild(this.recentItemsRoot);
    this.subwidgetsAdded.set(SpotlightTabNames.RECENT, this.recentItemsRoot);

    const adminActivity = this.recentUserActivity.filter((l) => l.description.toLowerCase().includes('published version'));
    this.adminRoot = adminActivity.length > 0 ? createSection(SpotlightTabNames.ADMIN, adminActivity, ui.icons.settings(() => {})) : null;
    if (this.adminRoot)
      root.appendChild(this.adminRoot);
    this.subwidgetsAdded.set(SpotlightTabNames.ADMIN, this.adminRoot);
    // }

    if (root.children.length === 0)
      root = this.getNewUserInfoColumns();

    const additionalFuncs = DG.Func.find({meta: {'isActivityWidget': 'true'}});
    for (const func of additionalFuncs) {
      const elements: HTMLElement[] = await func.apply();
      const list = ui.list(elements);
      list.classList.add('power-pack-activity-widget-subwidget-list-content');
      const listRoot = ui.divV([ui.h3(ui.span([ui.span([func.options['activityWidgetHeader'] ?? ''])]), 'power-pack-activity-widget-spotlight-column-header'),
        list], 'power-pack-activity-widget-spotlight-column');
      root.appendChild(listRoot);
    }
    const randomTip = this.tipsOfTheDay[Math.floor(Math.random() * this.tipsOfTheDay.length)];
    const tip = ui.divText(`💡 Tip of the day: ${randomTip}`, 'power-pack-activity-widget-spotlight-tip');
    root.appendChild(tip);

    setTimeout(() => this.cleanLists(), 500);
    console.timeEnd('ActivityDashboardWidget.buildSpotlightTab');
    return root;
  }

  getNewUserInfoColumns(): HTMLDivElement {
    const root = ui.divH([], 'power-pack-activity-widget-spotlight-root');
    const tutorialsApp = DG.Func.find({tags: ['app'], package: 'Tutorials', name: 'trackOverview'})[0];
    const demoApp = DG.Func.find({tags: ['app'], package: 'Tutorials', name: 'demoApp'})[0];
    if (!tutorialsApp || !demoApp)
      return root;
    const appHandler = DG.ObjectHandler.forEntity(demoApp);

    const gettingStartedList = ui.list([
      ui.link('Data Access', async () => await tutorialsApp.apply()),
      ui.link('Exploratory Data Analysis', async () => await tutorialsApp.apply()),
      ui.link('Cheminformatics', async () => await tutorialsApp.apply()),
      ui.link('Explore more Tutorials', async () => await tutorialsApp.apply()),
    ]);
    gettingStartedList.classList.add('power-pack-activity-widget-subwidget-list-content');
    const gettingStartedRoot = ui.divV([ui.h3(ui.span([ui.span(['Getting Started'])]), 'power-pack-activity-widget-spotlight-column-header'),
      gettingStartedList], 'power-pack-activity-widget-spotlight-column');

    const tryDemoAppsList = ui.list([
      ui.link('Scatter Plot', async () => await demoApp.apply({path: '/Visualization/General/Scatter-Plot'})),
      ui.link('Files', async () => await demoApp.apply({path: '/Data-Access/Files'})),
      ui.link('Matched Molecular Pairs', async () => await demoApp.apply({path: '/Cheminformatics/Matched-Molecular-Pairs'})),
      ui.link('Curve Fitting', async () => await demoApp.apply({path: '/Curves/Curve-Fitting'})),
    ]);
    tryDemoAppsList.classList.add('power-pack-activity-widget-subwidget-list-content');
    const tryDemoAppsRoot = ui.divV([ui.h3(ui.span([ui.span(['Try Demo Apps'])]), 'power-pack-activity-widget-spotlight-column-header'),
      tryDemoAppsList], 'power-pack-activity-widget-spotlight-column');

    let featuredApps: HTMLElement[] = [];
    const hitDesignApp = DG.Func.find({tags: ['app'], package: 'HitTriage', name: 'hitTriageApp'})[0];
    if (hitDesignApp)
      featuredApps.push((appHandler?.renderTooltip(DG.toDart(hitDesignApp)).firstChild as HTMLElement) ?? null);
    const admeticaApp = DG.Func.find({tags: ['app'], package: 'Admetica', name: 'admeticaApp'})[0];
    if (admeticaApp)
      featuredApps.push((appHandler?.renderTooltip(DG.toDart(admeticaApp)).firstChild as HTMLElement) ?? null);
    const diffStudioApp = DG.Func.find({tags: ['app'], package: 'DiffStudio', name: 'runDiffStudio'})[0];
    if (diffStudioApp)
      featuredApps.push((appHandler?.renderTooltip(DG.toDart(diffStudioApp)).firstChild as HTMLElement) ?? null);
    featuredApps = featuredApps.filter((app) => !!app);

    const featuredAppsList = ui.list([
      ...featuredApps,
      ui.link('Explore more Apps', () => setTimeout(() => grok.shell.addView(DG.View.createByType(DG.VIEW_TYPE.APPS)), 100)),
    ]);
    featuredAppsList.classList.add('power-pack-activity-widget-subwidget-list-content');
    const featuredAppsRoot = ui.divV([ui.h3(ui.span([ui.span(['Featured Apps'])]), 'power-pack-activity-widget-spotlight-column-header'),
      featuredAppsList], 'power-pack-activity-widget-spotlight-column');

    root.append(gettingStartedRoot, tryDemoAppsRoot, featuredAppsRoot);
    return root;
  }

  replaceMarkupSpans(container: HTMLElement, replacements: HTMLElement[], names: string[], time: dayjs.Dayjs): HTMLElement {
    const result = ui.span([], 'd4-markup');
    let replacementIndex = 0;
    for (let i = 0; i < container.childNodes.length; i++) {
      const node = container.childNodes[i];
      if (node instanceof HTMLSpanElement && node.textContent?.trim()?.startsWith('#{x')) {
        if (replacementIndex < replacements.length) {
          const span = ui.span([ui.span([
            ui.span([replacements[replacementIndex], ui.label(names[replacementIndex])], replacementIndex === 0 ? 'grok-markup-user' : ''),
          ], 'd4-link-label')]);
          result.appendChild(span);
          replacementIndex++;
        }
      } else {
        if (node instanceof Text && node.textContent === '.' && time)
          node.textContent = ` ${ui.time.timeSpan(time).textContent}`;
        result.appendChild(node);
      }
    }
    return result;
  }

  async getFavoritesTab(): Promise<HTMLElement> {
    console.time('ActivityDashboardWidget.buildFavoritesTab');
    const root = ui.div([]);
    const favorites = grok.shell.favorites;
    if (favorites.length === 0) {
      root.appendChild(ui.divText('No favorites found.'));
      return root;
    }

    const eventFetches = await Promise.all(favorites.map((entity) =>
      grok.dapi.log.where({entityId: entity.id,
        start: this.cutoffDate})
        .list()));
    this.favoritesEvents = (this.removeUnnecessaryEntities(eventFetches.flat()) as DG.LogEvent[])
      .sort((a, b) => b.eventTime?.diff(a.eventTime) ?? 0);

    if (this.favoritesEvents.length === 0) {
      root.appendChild(ui.divText('No recent activity on your favorites.'));
      return root;
    }
    this.favoritesListRoot = ui.list(this.favoritesEvents);
    root.appendChild(this.favoritesListRoot);
    this.cleanLists();
    console.timeEnd('ActivityDashboardWidget.buildFavoritesTab');
    return root;
  }

  async getNotificationsTab(): Promise<HTMLElement> {
    console.time('ActivityDashboardWidget.buildNotificationsTab');
    const root = ui.div([]);
    if (this.recentNotifications.length === 0) {
      root.appendChild(ui.divText('No recent notifications.'));
      return root;
    }

    this.notificationsListRoot = ui.list(this.recentNotifications);
    for (let i = 0; i < this.notificationsListRoot.children.length; i++) {
      const child = this.notificationsListRoot.children[i];
      const item = this.recentNotifications[i];
      if (item.createdAt) {
        const timestamp = ui.time.shortTimestamp(item.createdAt);
        timestamp.style.top = '5px';
        child.prepend(timestamp);
        const timeChild = child.querySelector('.d4-time');
        if (timeChild)
          timeChild.remove();
      }
    }
    root.appendChild(this.notificationsListRoot);
    this.cleanLists();
    console.timeEnd('ActivityDashboardWidget.buildNotificationsTab');
    return root;
  }

  async getActivityTab(): Promise<HTMLElement> {
    console.time('ActivityDashboardWidget.buildActivityTab');
    const recentUserActivityDataSource: DG.ActivityDataSource = grok.dapi.log.activity
      .where({userId: DG.User.current().id, start: this.cutoffDate});
    const recentActivity = await recentUserActivityDataSource.list();
    this.recentUserActivity = (this.removeUnnecessaryEntities(recentActivity
      .filter((a) => a.eventTime?.isAfter(this.cutoffDate))) as DG.LogEvent[])
      .sort((a, b) => b.eventTime?.diff(a.eventTime) ?? 0);

    const root = ui.div([]);
    if (this.recentUserActivity.length === 0) {
      root.appendChild(ui.divText('No recent user activity.'));
      return root;
    }
    this.activityListRoot = ui.list(this.recentUserActivity);
    root.appendChild(this.activityListRoot);
    this.cleanLists();
    console.timeEnd('ActivityDashboardWidget.buildActivityTab');
    return root;
  }

  cleanLists(): void {
    if (this.favoritesListRoot && this.favoritesListRoot.children.length > 0)
      this.cleanList(this.favoritesListRoot);
    if (this.notificationsListRoot && this.notificationsListRoot.children.length > 0)
      this.cleanList(this.notificationsListRoot, false);
    if (this.activityListRoot && this.activityListRoot.children.length > 0)
      this.cleanList(this.activityListRoot);
    if (this.actionRequiredRoot && this.actionRequiredRoot.children.length > 0 && this.actionRequiredRoot.lastElementChild) {
      const list = this.actionRequiredRoot.lastElementChild;
      for (let i = 0; i < list.children.length; i++) {
        const listChild = list.children[i];
        const text = listChild.textContent?.toLowerCase();
        if ((text?.includes('you were assigned report #') || text?.includes('you were assigned  report #')) && listChild.querySelector('.d4-markup label')) {
          listChild.querySelector('.d4-markup label')!.textContent = `${this.reportAmount} ${this.reportAmount && this.reportAmount > 1 ? 'reports' : 'report'}`;
          const reportsApp = DG.Func.find({package: 'UsageAnalysis', name: 'reportsApp'})[0];
          if (reportsApp) {
            const handler = DG.ObjectHandler.forEntity(reportsApp);
            const reportsAppLabelSpan = handler?.renderTooltip(DG.toDart(reportsApp)).firstChild;
            if (reportsAppLabelSpan && listChild.lastElementChild)
              listChild.lastElementChild.appendChild(ui.span(['. Open ', reportsAppLabelSpan, ' to see more.']));
          }
        }
      }
    }
  }

  cleanList(list: HTMLElement, aggregateUnique: boolean = true): HTMLElement {
    const uniqueEvents: Map<string, number> = new Map<string, number>();
    for (let i = 0; i < list.children.length; i++) {
      const child = list.children[i];
      const text = child.textContent;
      if (aggregateUnique && text) {
        if (uniqueEvents.has(text)) {
          uniqueEvents.set(text, uniqueEvents.get(text)! + 1);
          child.remove();
          continue;
        }
        else
          uniqueEvents.set(text, 1);
      }
      const div = (child as HTMLElement).querySelector('div');
      if (div) {
        if (div.parentElement && div.children.length > 0)
          div.parentElement.insertBefore(ui.span(Array.from(div.children)), div);
        div.remove();
      }
      child.classList.add('power-pack-activity-widget-row');
    }
    for (let i = 0; i < list.children.length; i++) {
      const child = list.children[i];
      const text = child.textContent;
      if (!text)
        continue;
      if (uniqueEvents.has(text) && !text.endsWith(' times') && uniqueEvents.get(text)! > 1)
        child.querySelector('span.d4-markup')?.appendChild(ui.span([` ${uniqueEvents.get(text)} times`]));
    }
    return list;
  }

  removeUnnecessaryEntities(list: Array<DG.LogEvent | DG.UserNotification>): Array<DG.LogEvent | DG.UserNotification> {
    return list.filter((item) => {
      const text = item instanceof DG.LogEvent ? item.description?.toLowerCase() : item.text?.toLowerCase();
      return !this.keywordsToIgnore.some((keyword) => text?.includes(keyword));
    });
  }
}
