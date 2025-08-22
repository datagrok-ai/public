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
  static RECORDS_AMOUNT_PER_PAGE = 50;

  keywordsToIgnore: string[] = ['"data"', '"data w/', '"fitted data"', '"reduced data"', '"reduceddata"', '"sizingparams"', '"primaryfilter', '"trimmed data"', '"input data from imported file"', '"dataanalysisdf', '"clean data"'/*, 'Ran <span>#{x.'*/];
  currentDate: dayjs.Dayjs = dayjs();

  favoritesEvents: DG.LogEvent[] = [];

  recentNotifications: DG.UserNotification[] = [];
  recentUserActivity: DG.LogEvent[] = [];

  sharedNotifications: DG.UserNotification[] = [];
  sharedUsers: DG.User[] = [];
  sharedWithMe: DG.Entity[] = [];
  recentEntities: DG.Entity[] = [];
  recentEntityTimes: dayjs.Dayjs[] = [];

  spotlightRoot?: HTMLElement;
  favoritesListRoot?: HTMLElement;
  notificationsListRoot?: HTMLElement;
  activityListRoot?: HTMLElement;

  actionRequiredRoot?: HTMLDivElement | null = null;
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

    const tabControl = ui.tabControl(tabs, true);
    tabControl.onTabChanged.subscribe((_) => this.cleanLists());
    tabControl.root.style.height = '100%';
    tabControl.root.style.width = '100%';

    this.root.appendChild(tabControl.root);
    setTimeout(() => {
      const header = this.root.parentElement?.parentElement?.querySelector('.d4-dialog-header');
      if (header) {
        const icons = Array.from(header.getElementsByClassName('grok-icon'));
        this.root.appendChild(icons[icons.length - 1]);
        header.remove();
      }
    }, 500);
  }

  async initSpotlightData(): Promise<void> {
    console.time('ActivityDashboardWidget.notifications');
    const notificationsDataSource: DG.NotificationsDataSource = grok.dapi.users.notifications.forCurrentUser()
      .by(ActivityDashboardWidget.RECORDS_AMOUNT_PER_PAGE) as DG.NotificationsDataSource;
    const notifications = await notificationsDataSource.list();
    this.recentNotifications = notifications
      .filter((n) => !n.isRead || n.createdAt?.isAfter(this.currentDate.subtract(ActivityDashboardWidget.RECENT_TIME_DAYS, 'day')))
      .sort((a, b) => b.createdAt?.diff(a.createdAt) ?? 0);
    this.sharedNotifications = this.recentNotifications.filter((n) => n.text.startsWith('<span>#') &&
      n.text.includes('</span> shared') && n.text.includes(': <span>'));
    const sharedUserIds: string[] = [];
    const sharedEntityIds: string[] = [];
    this.sharedNotifications.forEach((n) => {
      const matches = [...n.text.matchAll(/<span>#{x\.([a-f0-9-]+)\.".*?"}<\/span>/g)];
      if (matches.length > 1) {
        sharedUserIds.push(matches[0][1]);
        sharedEntityIds.push(matches[1][1]);
      }
    });
    console.timeEnd('ActivityDashboardWidget.notifications');
    console.time('ActivityDashboardWidget.getEntitiesByIds');
    const entities = await grok.dapi.getEntities([...sharedUserIds, ...sharedEntityIds]);
    this.sharedUsers = entities.filter((ent) => ent instanceof DG.User) as DG.User[];
    this.sharedWithMe = entities.filter((ent) => !(ent instanceof DG.User));
    console.timeEnd('ActivityDashboardWidget.getEntitiesByIds');

    console.time('ActivityDashboardWidget.mostRecentEntities');
    const mostRecentEntitiesDf = await queries.mostRecentEntities(DG.User.current().id);
    console.timeEnd('ActivityDashboardWidget.mostRecentEntities');
    const recentEntityIdCol = mostRecentEntitiesDf.col('id');
    const lastEventTimeCol = mostRecentEntitiesDf.col('last_event_time');
    if (recentEntityIdCol && lastEventTimeCol) {
      const recentEntityIds = new Array(recentEntityIdCol.length);
      console.time('ActivityDashboardWidget.mostRecentEntities.getAllEntities');
      for (let i = 0; i < recentEntityIdCol.length; i++)
        recentEntityIds[i] = recentEntityIdCol.get(i);
      const recentEntsNotFiltered = await grok.dapi.getEntities(recentEntityIds);
      for (let i = 0; i < recentEntsNotFiltered.length; i++) {
        const ent = recentEntsNotFiltered[i];
        if (!(ent instanceof DG.FuncCall || ent instanceof DG.Group || ent instanceof DG.User || ent instanceof DG.Package ||
          ent instanceof DG.UserReport || ent instanceof DG.TableInfo || (ent instanceof DG.Func && !(ent instanceof DG.Script ||
          ent instanceof DG.DataQuery || ent instanceof DG.DataJob)) || ent instanceof DG.ViewInfo || ent == null ||
          //@ts-ignore
          (ent instanceof DG.Project && (!ent.isDashboard || ent.isPackage)) || (ent.hasOwnProperty('npmScope') && ent['npmScope'] == 'datagrok'))) {
          this.recentEntities.push(ent);
          this.recentEntityTimes.push(lastEventTimeCol.get(i));
        }
      }
      console.timeEnd('ActivityDashboardWidget.mostRecentEntities.getAllEntities');
    }
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
      const ITEMS_LENGTH = 10;
      const itemsToShow = usedItems.slice(0, ITEMS_LENGTH);
      const list = ui.list(itemsToShow);
      list.classList.add('power-pack-activity-widget-subwidget-list-content');
      const listChildren = Array.from(list.children);
      for (let i = 0; i < itemsToShow.length; i++) {
        const item = itemsToShow[i];
        const listChild = listChildren[i];
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
      return ui.divV([ui.h3(ui.span([icon, ui.span([` ${title}`])]), 'power-pack-activity-widget-spotlight-column-header'),
        list], title === SpotlightTabNames.ACTION_REQUIRED ? 'power-pack-activity-widget-spotlight-column-action-required' : 'power-pack-activity-widget-spotlight-column');
    };

    const root = ui.divH([], 'power-pack-activity-widget-spotlight-root');

    const actionRequired = this.recentNotifications.filter((n) => {
      const text = n.text.toLowerCase();
      return text.includes('you were assigned') || text.includes('requested a membership');
    });
    this.actionRequiredRoot = actionRequired.length > 0 ? createSection(SpotlightTabNames.ACTION_REQUIRED, actionRequired, ui.iconFA('exclamation-circle')) : null;
    if (this.actionRequiredRoot) {
      if (actionRequired.some((n) => n.text.includes('requested a membership')))
        this.actionRequiredRoot.style.maxWidth = '200px';
      root.appendChild(this.actionRequiredRoot);
    }

    const sharedWithMeRoot = this.sharedWithMe.length > 0 ? createSection(SpotlightTabNames.SHARED_WITH_ME, this.sharedWithMe, ui.iconFA('inbox')) : null;
    if (sharedWithMeRoot)
      root.appendChild(sharedWithMeRoot);

    const spacesRoot = null;
    if (spacesRoot)
      root.appendChild(spacesRoot);

    const recentItemsRoot = this.recentEntities.length > 0 ? createSection(SpotlightTabNames.RECENT, this.recentEntities, ui.iconFA('history')) : null;
    if (recentItemsRoot)
      root.appendChild(recentItemsRoot);

    const adminActivity = this.recentUserActivity.filter((l) => l.description.toLowerCase().includes('published version'));
    const adminRoot = adminActivity.length > 0 ? createSection(SpotlightTabNames.ADMIN, adminActivity, ui.icons.settings(() => {})) : null;
    if (adminRoot)
      root.appendChild(adminRoot);
    this.cleanLists();

    console.timeEnd('ActivityDashboardWidget.buildSpotlightTab');

    return root;
  }

  replaceMarkupSpans(container: HTMLElement, replacements: HTMLElement[], names: string[], time: dayjs.Dayjs): HTMLElement {
    const result = ui.span([], 'd4-markup');
    let replacementIndex = 0;
    const children = Array.from(container.childNodes);
    for (let i = 0; i < children.length; i++) {
      const node = children[i];
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
        start: this.currentDate.subtract(ActivityDashboardWidget.RECENT_TIME_DAYS, 'days')}) // where entityId in
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
    const childList = Array.from(this.notificationsListRoot.children);
    for (let i = 0; i < childList.length; i++) {
      const child = childList[i];
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
      .where({userId: DG.User.current().id, start: this.currentDate.subtract(ActivityDashboardWidget.RECENT_TIME_DAYS, 'day')});
    const recentActivity = await recentUserActivityDataSource.list();
    this.recentUserActivity = (this.removeUnnecessaryEntities(recentActivity
      .filter((a) => a.eventTime?.isAfter(this.currentDate.subtract(ActivityDashboardWidget.RECENT_TIME_DAYS, 'day')))) as DG.LogEvent[])
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
        if (listChild.textContent?.toLowerCase()?.includes('you were assigned  report #'))
          listChild.querySelector('.d4-markup label')!.textContent = `${this.reportAmount} reports`;
      }
    }
  }

  cleanList(list: HTMLElement, aggregateUnique: boolean = true): HTMLElement {
    const uniqueEvents: Map<string, number> = new Map<string, number>();
    const oldChildren = Array.from(list.children);
    for (const child of oldChildren) {
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
    const newChildren = Array.from(list.children);
    for (const child of newChildren) {
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
