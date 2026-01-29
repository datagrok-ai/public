import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';
import {queries} from '../package-api';
import {LearningWidget} from './learning-widget';


enum SpotlightTabNames {
  ACTION_REQUIRED = 'Action required',
  SHARED_WITH_ME = 'Shared with me',
  RECENT = 'Recent',
}

export class ActivityDashboardWidget extends DG.Widget {
  get type(): string {
    return 'ActivityDashboardWidget';
  }

  static RECENT_TIME_DAYS = 2;
  static SPOTLIGHT_ITEMS_LENGTH = 8;

  static sharedEntityRegex: RegExp = /<span>#{x\.([a-f0-9-]+)\.".*?"}<\/span>/g;
  keywordsToIgnore: string[] = ['"data"', '"data w/', '"fitted data"', '"reduced data"', '"reduceddata"', '"sizingparams"', '"primaryfilter', '"trimmed data"', '"input data from imported file"', '"dataanalysisdf', '"clean data"'/*, 'Ran <span>#{x.'*/];
  tipsOfTheDay: string[] = [
    'To sort the data, double-click a column header.',
    'To load a CSV or Excel file as a table, drag and drop it into Datagrok.',
    'To see statistics of a column, hover over its name.',
    'To select multiple columns, hold Shift and click several column headers.',
    'To search across all values in a table, use Ctrl+F.',
    'To resize columns, drag the separator in the header.',
    'To save and restore your favorite viewer arrangements, use Layouts.',
    'To create calculated columns with formulas, use “Add New Column”.',
    'To see context actions for a column, right-click it.',
    'To quickly explore relationships between numeric columns, try the “Correlation Plot”.',
  ];
  availableDemosOfTheDay: DG.Func[] = DG.Func.find({meta: {'demoPath': null}});
  tutorialsOfTheDay: string[] = ['Grid Customization', 'Viewers', 'Scatter Plot', 'Embedded Viewers', 'Filters', 'Dashboards',
    'Multivariate Analysis', 'Scripting', 'R-Groups Analysis', 'Activity Cliffs', 'Similarity and Diversity Search',
    'Substructure Search and Filtering', 'Data Connectors', 'Data Aggregation', 'Calculated Columns', 'Differential equations',
    'Sensitivity analysis', 'Parameter optimization'];
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
  subwidgetsAdded: Map<string, HTMLDivElement | null | undefined> = new Map<string, HTMLDivElement>();
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
      'Learn': () => new LearningWidget().root,
    };

    this.tabControl = ui.tabControl(tabs, true);
    this.subs.push(this.tabControl.onTabChanged.subscribe((tabPane: DG.TabPane) => {
      this.cleanLists();
      tabPane.name === 'Learn' ? tabPane.content.parentElement?.classList.add('power-pack-overflow-hidden') :
        tabPane.content.parentElement?.classList.remove('power-pack-overflow-hidden');
    }));
    this.tabControl.root.style.height = '90%';
    this.tabControl.root.style.width = '100%';

    this.root.appendChild(this.tabControl.root);
    this.root.appendChild(await this.createRandomizedTipOfTheDay());
    // TODO: uncomment when JS API 1.27/1.26.5 is released
    // this.subs.push((grok.events.onCurrentObjectChanged.subscribe((args) => {
    //   const trigger = args.args.trigger;
    //   if (this.root.contains(trigger))
    //     grok.shell.windows.context.visible = true;
    // })));
  }

  async getDemosOfTheDay(): Promise<string[]> {
    const demoAppHierarchyFunc = DG.Func.find({name: 'getDemoAppHierarchy'})[0];
    if (demoAppHierarchyFunc == null)
      return [];
    const demoAppHierarchy = await demoAppHierarchyFunc.apply();
    if (demoAppHierarchy == null || demoAppHierarchy.length === 0)
      return [];
    const parsedHierarchy = JSON.parse(demoAppHierarchy);
    const flattenDemoHierarchy = (node: any, path: string[] = []) => {
      const currentPath = [...path, node.name];
      if (!node.children || node.children.length === 0)
        return [currentPath.join(' | ')]; // full path "Category | Subcategory | Demo"
      return node.children.flatMap((child: any) => flattenDemoHierarchy(child, currentPath));
    };
    return parsedHierarchy.children.flatMap((child: any) => flattenDemoHierarchy(child));
  }

  async createRandomizedTipOfTheDay(): Promise<HTMLElement> {
    const seededRandom = (seed: string) => {
      let x = 0;
      for (let i = 0; i < seed.length; i++) {
        x += seed.charCodeAt(i);
        x = (x * 9301 + 49297) % 233280;
      }
      return x / 233280;
    };
    const getISOWeek = (date: Date) => {
      const tmpDate = new Date(date.getTime());
      tmpDate.setHours(0, 0, 0, 0);
      tmpDate.setDate(tmpDate.getDate() + 3 - ((tmpDate.getDay() + 6) % 7)); // Thursday of this week
      const week1 = new Date(tmpDate.getFullYear(), 0, 4); // first Thursday of the year
      return 1 + Math.round(((tmpDate.getTime() - week1.getTime()) / 86400000 - 3 + ((week1.getDay() + 6) % 7)) / 7);
    };
    const today = new Date();
    const year = today.getFullYear();
    const weekNumber = getISOWeek(today);
    const dayNames = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'];
    const dayName = [0, 6].includes(today.getDay()) ? 'Friday' : dayNames[today.getDay() - 1];
    const weekSeed = `${year}-W${weekNumber}-${dayName}`;

    const dayTipTypeMap: Record<string, 'tip' | 'demo' | 'tutorial'> = {
      Monday: 'demo',
      Tuesday: 'tip',
      Wednesday: 'tutorial',
      Thursday: 'tip',
      Friday: 'demo',
      Saturday: 'demo', // repeat Friday
      Sunday: 'demo', // repeat Friday
    };
    const tipType = dayTipTypeMap[dayName];
    const todayItemListFromType = tipType === 'tip' ? this.tipsOfTheDay : tipType === 'tutorial' ? this.tutorialsOfTheDay : [];
    const randomizedTips = [...this.tipsOfTheDay, ...this.availableDemosOfTheDay, ...this.tutorialsOfTheDay];
    let randomTip: DG.Func | string;
    if (tipType === 'demo') {
      let demosOfTheDay: (string | DG.Func)[] = await this.getDemosOfTheDay();
      if (demosOfTheDay.length === 0)
        demosOfTheDay = this.availableDemosOfTheDay;
      randomTip = demosOfTheDay.length > 0 ? demosOfTheDay[Math.floor(seededRandom(weekSeed) * demosOfTheDay.length)] : '';
      if (!(randomTip instanceof DG.Func)) {
        const demoFunc = DG.Func.find({meta: {'demoPath': randomTip}})[0];
        randomTip = demoFunc ? demoFunc : this.availableDemosOfTheDay[Math.floor(seededRandom(weekSeed) * this.availableDemosOfTheDay.length)];
      }
    }
    else
      randomTip = todayItemListFromType.length > 0 ? todayItemListFromType[Math.floor(seededRandom(weekSeed) * todayItemListFromType.length)] : '';

    let tipIdx = 0;
    const demoApp = DG.Func.find({meta: {role: DG.FUNC_TYPES.APP}, package: 'Tutorials', name: 'demoApp'})[0];
    const tutorialsApp = DG.Func.find({meta: {role: DG.FUNC_TYPES.APP}, package: 'Tutorials', name: 'trackOverview'})[0];

    const createTip = (newRandomTip: DG.Func | string) => {
      const tip = ui.divText('', 'power-pack-activity-widget-spotlight-tip');
      const tipText = ui.span([`💡 ${newRandomTip instanceof DG.Func ? 'Demo' : this.tutorialsOfTheDay.includes(newRandomTip) ?
        'Tutorial' : 'Tip'} of the day: ${newRandomTip instanceof DG.Func || this.tutorialsOfTheDay.includes(newRandomTip) ? '' : newRandomTip}`]);
      tip.appendChild(tipText);
      let link: HTMLAnchorElement;
      if (newRandomTip instanceof DG.Func && demoApp) {
        const path = newRandomTip.options[DG.FUNC_OPTIONS.DEMO_PATH] as string;
        const pathArray = path.split('|').map((s) => s.trim());
        const actualPath = pathArray.map((s) => s.replaceAll(' ', '-')).join('/');
        link = ui.link(pathArray[pathArray.length - 1], async () => await demoApp.apply({path: `/${actualPath}`}), newRandomTip.description);
        tip.appendChild(link);
      }
      else if (!(newRandomTip instanceof DG.Func) && this.tutorialsOfTheDay.includes(newRandomTip) && tutorialsApp) {
        link = ui.link(newRandomTip, async () => await tutorialsApp.apply());
        tip.appendChild(link);
      }
      const updateTip = () => {
        tipIdx++;
        const newRandomTip = randomizedTips[Math.floor(Math.random() * randomizedTips.length)];
        const newTip = createTip(newRandomTip);
        tip.replaceWith(newTip);
      };

      const nextTipIcon = ui.iconFA('chevron-right', updateTip, 'Next tip');
      nextTipIcon.classList.add('tip-next');
      tip.prepend(nextTipIcon);
      return tip;
    };
    return createTip(randomTip);
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

    const recentEntityIdCol = mostRecentEntitiesDf?.col('id')!;
    const lastEventTimeCol = mostRecentEntitiesDf?.col('last_event_time')!;
    const recentEntityIds = new Array(recentEntityIdCol?.length ?? 0);
    if (recentEntityIdCol)
      for (let i = 0; i < recentEntityIdCol.length; i++)
        recentEntityIds[i] = recentEntityIdCol.get(i);

    const uniqueIds = Array.from(new Set<string>([...sharedUserIds, ...sharedEntityIds, ...recentEntityIds]));
    console.time('ActivityDashboardWidget.getEntitiesByIds');
    const allEntities = await grok.dapi.getEntities(uniqueIds);
    console.timeEnd('ActivityDashboardWidget.getEntitiesByIds');

    const byIdMap = new Map(allEntities.filter((e) => e != null).map((e) => [e.id, e]));

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
        (ent instanceof DG.Project && (!ent.isDashboard || ent.isPackage)) || (ent.hasOwnProperty('npmScope') && ent['npmScope'] == 'datagrok')) && lastEventTimeCol && lastEventTimeCol.get(i)) {
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
        if (item instanceof DG.Notebook || item instanceof DG.Model || item instanceof DG.ViewLayout) {
          const parentElem = listChild.firstChild;
          const iconToPaste = item instanceof DG.Notebook ? ui.iconImage('notebook', '/images/entities/jupyter.png') :
            item instanceof DG.Model ? ui.iconSvg('model') : ui.iconSvg('view-layout');
          if (parentElem && parentElem.firstChild)
            parentElem.insertBefore(iconToPaste, parentElem.firstChild);
        } else if (item instanceof DG.UserNotification) {
          if (item.createdAt) {
            const timestamp = ui.time.shortTimestamp(item.createdAt);
            timestamp.style.top = '3px';
            listChild.prepend(timestamp);
            const timeChild = listChild.querySelector('.d4-time');
            if (timeChild)
              timeChild.remove();
          }
        }

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
          icon.style.top = '3px';
          listChild.prepend(icon);
        } else if (title === SpotlightTabNames.RECENT) {
          const timestamp = ui.time.shortTimestamp(this.recentEntityTimes[i]);
          timestamp.style.top = '3px';
          listChild.prepend(timestamp);
        }
      }
      this.subwidgetsAmount++;
      return ui.divV([ui.h3(ui.span([icon, ui.span([` ${title}`])]), 'power-pack-activity-widget-spotlight-column-header'),
        list], title === SpotlightTabNames.ACTION_REQUIRED ? 'power-pack-activity-widget-spotlight-column-action-required' : 'power-pack-activity-widget-spotlight-column');
    };

    let root = ui.divH([], 'power-pack-activity-widget-spotlight-root');

    if (!(DG.User.current().joined > dayjs().subtract(5, 'day'))) {
      const actionRequired = this.recentNotifications.filter((n) => {
        const text = n.text.toLowerCase();
        return text.includes('you were assigned') || text.includes('requested a membership');
      });
      this.actionRequiredRoot = actionRequired.length > 0 ? createSection(SpotlightTabNames.ACTION_REQUIRED, actionRequired, ui.iconFA('exclamation-circle')) : null;
      if (this.actionRequiredRoot)
        root.appendChild(this.actionRequiredRoot!);
      this.subwidgetsAdded.set(SpotlightTabNames.ACTION_REQUIRED, this.actionRequiredRoot);

      this.sharedWithMeRoot = this.sharedWithMe.length > 0 ? createSection(SpotlightTabNames.SHARED_WITH_ME, this.sharedWithMe, ui.iconFA('inbox')) : null;
      if (this.sharedWithMeRoot)
        root.appendChild(this.sharedWithMeRoot!);
      this.subwidgetsAdded.set(SpotlightTabNames.SHARED_WITH_ME, this.sharedWithMeRoot);

      this.spacesRoot = null;
      if (this.spacesRoot)
        root.appendChild(this.spacesRoot!);
      // this.subwidgetsAdded.set(SpotlightTabNames.SPACES, this.spacesRoot);

      this.recentItemsRoot = this.recentEntities.length > 0 ? createSection(SpotlightTabNames.RECENT, this.recentEntities, ui.iconFA('history')) : null;
      if (this.recentItemsRoot)
        root.appendChild(this.recentItemsRoot!);
      this.subwidgetsAdded.set(SpotlightTabNames.RECENT, this.recentItemsRoot);
    }

    if (root.children.length === 0)
      root = await this.getNewUserInfoColumns();

    const additionalFuncs = DG.Func.find({meta: {'activityWidgetHeader': null}});
    for (const func of additionalFuncs) {
      const subWidget: DG.Widget = await func.apply();
      subWidget.root.classList.add('power-pack-activity-widget-subwidget-list-content');
      const rootToAppend = ui.divV([ui.h3(ui.span([ui.span([func.options['activityWidgetHeader'] ?? ''])]), 'power-pack-activity-widget-spotlight-column-header'),
        subWidget.root], 'power-pack-activity-widget-spotlight-column');
      root.appendChild(rootToAppend);
    }

    setTimeout(() => this.cleanLists(), 500);
    console.timeEnd('ActivityDashboardWidget.buildSpotlightTab');
    return root;
  }

  async getNewUserInfoColumns(): Promise<HTMLDivElement> {
    const root = ui.divH([], 'power-pack-activity-widget-spotlight-root');
    const tutorialsApp = DG.Func.find({meta: {role: DG.FUNC_TYPES.APP}, package: 'Tutorials', name: 'trackOverview'})[0];
    const demoApp = DG.Func.find({meta: {role: DG.FUNC_TYPES.APP}, package: 'Tutorials', name: 'demoApp'})[0];
    if (!tutorialsApp || !demoApp)
      return root;
    const appHandler = DG.ObjectHandler.forEntity(demoApp);

    const createLinkWithIcon = (text: string, target: Function, appFunc: DG.Func) => {
      const link = ui.link('', target);
      const appElement = appHandler?.renderTooltip(DG.toDart(appFunc)).firstChild as HTMLElement;
      if (appElement) {
        const appElementLabel = appElement.querySelector('label');
        if (appElementLabel)
          appElementLabel.textContent = text;
        link.appendChild(appElement);
      }
      return link;
    };

    const getDemoDashboardName = (demo: DG.Func) => {
      const path = demo.options[DG.FUNC_OPTIONS.DEMO_PATH] as string;
      const pathArray = path.split('|').map((s) => s.trim());
      return pathArray[pathArray.length - 1];
    };

    const sortDemoDashboards = (items: DG.Func[], priorities: string[]): DG.Entity[] => {
      const getPriority = (item: DG.Func) => {
        for (let i = 0; i < priorities.length; i++)
          if (getDemoDashboardName(item).toLowerCase().includes(priorities[i].toLowerCase()))
            return i;
        return priorities.length;
      };
      return items.map((item, index) => ({item, index})).sort((a, b) => {
        const pa = getPriority(a.item);
        const pb = getPriority(b.item);
        if (pa !== pb)
          return pa - pb;
        return a.index - b.index;
      }).map((obj) => obj.item);
    };

    const gettingStartedList = ui.list([
      createLinkWithIcon('Data Access', async () => await tutorialsApp.apply(), tutorialsApp),
      createLinkWithIcon('Data Transformation', async () => await tutorialsApp.apply(), tutorialsApp),
      createLinkWithIcon('Exploratory Data Analysis', async () => await tutorialsApp.apply(), tutorialsApp),
      createLinkWithIcon('Cheminformatics', async () => await tutorialsApp.apply(), tutorialsApp),
      createLinkWithIcon('Explore more Tutorials', async () => await tutorialsApp.apply(), tutorialsApp),
    ]);
    gettingStartedList.classList.add('power-pack-activity-widget-subwidget-list-content', 'power-pack-activity-widget-getting-started-subwidget-list-content');
    const gettingStartedRoot = ui.divV([ui.h3(ui.span([ui.iconFA('rocket'), ui.span([' Interactive Tutorials'])]), 'power-pack-activity-widget-spotlight-column-header'),
      gettingStartedList], 'power-pack-activity-widget-spotlight-getting-started-column');

    const tryDemoAppsList = ui.list([
      createLinkWithIcon('Scatter Plot', async () => await demoApp.apply({path: '/Visualization/General/Scatter-Plot'}), demoApp),
      createLinkWithIcon('Files', async () => await demoApp.apply({path: '/Data-Access/Files'}), demoApp),
      createLinkWithIcon('Matched Molecular Pairs', async () => await demoApp.apply({path: '/Cheminformatics/Matched-Molecular-Pairs'}), demoApp),
      createLinkWithIcon('Trellis Plot', async () => await demoApp.apply({path: '/Visualization/Data-Separation/Trellis-Plot'}), demoApp),
      createLinkWithIcon('Curve Fitting', async () => await demoApp.apply({path: '/Curves/Curve-Fitting'}), demoApp),
    ]);
    tryDemoAppsList.classList.add('power-pack-activity-widget-subwidget-list-content', 'power-pack-activity-widget-getting-started-subwidget-list-content');
    const tryDemoAppsRoot = ui.divV([ui.h3(ui.span([ui.iconFA('play-circle'), ui.span([' Demo Apps'])]), 'power-pack-activity-widget-spotlight-column-header'),
      tryDemoAppsList], 'power-pack-activity-widget-spotlight-getting-started-column');
    root.append(gettingStartedRoot, tryDemoAppsRoot);

    const dashboardsToShow: HTMLElement[] = [];
    let medChemDashboard: DG.Project | undefined | null = null;
    const medChemDashboardId = grok.userSettings.getValue('activity-dashboard-storage', 'medChemDashboardId');
    if (!medChemDashboardId) {
      const demoProjects = await grok.dapi.projects.filter('#demo').list();
      medChemDashboard = demoProjects.find((p) => p.name.toLowerCase().includes('medchem'));
      if (medChemDashboard)
        grok.userSettings.put('activity-dashboard-storage', {'medChemDashboardId': medChemDashboard.id});
    }
    else
      medChemDashboard = await grok.dapi.projects.find(medChemDashboardId);

    const otherDemoDashboards = DG.Func.find({meta: {'isDemoDashboard': 'true'}});
    const sortedDemoDashboards = sortDemoDashboards(otherDemoDashboards, ['peptide sar', 'sequence space', 'r-group analysis', 'molecule activity cliffs']);
    const demoDashboards: DG.Entity[] = [...(medChemDashboard ? [medChemDashboard] : []), ...sortedDemoDashboards];
    for (let i = 0; i < Math.min(5, demoDashboards.length); i++) {
      const demo = demoDashboards[i] instanceof DG.Project ? demoDashboards[i] as DG.Project : demoDashboards[i] as DG.Func;
      const link = demo instanceof DG.Project ? ui.link('', () => demo.open()) : ui.link('', async () => await demo.apply());
      link.append(ui.iconSvg('project'), ui.span([demo instanceof DG.Project ? demo.friendlyName : getDemoDashboardName(demo)]));
      dashboardsToShow.push(ui.span([link], 'd4-link-label'));
    }
    if (dashboardsToShow.length > 0) {
      const dashboardList = ui.list(dashboardsToShow);
      dashboardList.classList.add('power-pack-activity-widget-subwidget-list-content', 'power-pack-activity-widget-getting-started-subwidget-list-content');
      const dashboardsRoot = ui.divV([ui.h3(ui.span([ui.iconFA('table'), ui.span([' Demo Datasets'])]), 'power-pack-activity-widget-spotlight-column-header'),
        dashboardList], 'power-pack-activity-widget-spotlight-getting-started-column');
      root.append(dashboardsRoot);
    }
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
    const root = ui.div([], 'power-pack-activity-widget-favorites-tab');
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
    const root = ui.div([], 'power-pack-activity-widget-notifications-tab');
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

    const root = ui.div([], 'power-pack-activity-widget-activity-tab');
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
