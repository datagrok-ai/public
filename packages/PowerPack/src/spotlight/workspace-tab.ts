import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';
import {getMyGroupFavorites, unpin} from './group-favorites';
import type {SpotlightWidget} from './spotlight-widget';


/** Returns the platform icon for an entity via its registered EntityMeta handler. */
function entityIcon(entity: DG.Entity): HTMLElement {
  return DG.ObjectHandler.forEntity(entity)?.renderIcon(entity.dart) ?? ui.iconFA('file');
}

/** Opens an entity in the appropriate way. */
function openEntity(entity: DG.Entity): void {
  if (entity instanceof DG.Project)
    entity.open();
  else if (entity instanceof DG.Func)
    entity.apply();
  else
    grok.shell.o = entity;
}

/** Creates a single card tile for an entity. */
function createCard(entity: DG.Entity, root: HTMLElement, group?: DG.Group): HTMLElement {
  const icon = entityIcon(entity);
  const iconWrap = ui.div([icon], 'pp-workspace-card-icon');

  const label = ui.div([entity.friendlyName], 'pp-workspace-card-label');
  const card = ui.divV([iconWrap, label], 'pp-workspace-card');
  card.addEventListener('click', () => openEntity(entity));
  ui.tooltip.bind(card, entity.friendlyName);
  ui.bind(entity, card, {contextMenu: true});

  if (group) {
    const removeBtn = ui.icons.close(async () => {
      await unpin(entity, group);
      const section = card.closest('.pp-workspace-section')!;
      card.remove();
      if (section.querySelectorAll('.pp-workspace-card').length === 0) {
        section.remove();
        if (root.querySelectorAll('.pp-workspace-section').length === 0) {
          root.appendChild(ui.divText(
            'Your workspace is empty. Ask your admin to pin items for your group, or explore on your own.',
            'pp-workspace-empty-hint',
          ));
        }
      }
    }, `Unpin from ${group.friendlyName}`);
    removeBtn.classList.add('pp-workspace-card-remove');
    removeBtn.addEventListener('click', (e: MouseEvent) => e.stopPropagation());
    card.appendChild(removeBtn);
  }

  return card;
}

/** Creates a section with header and a row of cards. */
function createSection(title: string, iconName: string, entities: DG.Entity[], root: HTMLElement, group?: DG.Group): HTMLElement {
  const headerIcon = ui.iconFA(iconName);
  const header = ui.h3(
    ui.span([headerIcon, ui.span([` ${title}`])]),
    'pp-workspace-section-header',
  );
  const cards = entities.map((e) => createCard(e, root, group));
  const row = ui.divH(cards, 'pp-workspace-cards-row');
  return ui.divV([header, row], 'pp-workspace-section');
}

/** Creates a recent-items list section. */
function createRecentSection(
  entities: DG.Entity[],
  times: (dayjs.Dayjs | null)[],
): HTMLElement {
  const headerIcon = ui.iconFA('history');
  const header = ui.h3(
    ui.span([headerIcon, ui.span([' Continue where you left off'])]),
    'pp-workspace-section-header',
  );

  const items: HTMLElement[] = [];
  for (let i = 0; i < entities.length; i++) {
    const entity = entities[i];
    const icon = entityIcon(entity);
    icon.style.width = '16px';
    icon.style.height = '16px';

    const link = ui.link(entity.friendlyName, () => openEntity(entity));
    const row = ui.divH([icon, link], 'pp-workspace-recent-item');

    const time = times[i];
    if (time) {
      const timestamp = ui.time.shortTimestamp(time);
      timestamp.classList.add('pp-workspace-recent-time');
      row.appendChild(timestamp);
    }

    ui.bind(entity, row, {contextMenu: true});
    items.push(row);
  }

  const list = ui.divV(items, 'pp-workspace-recent-list');
  return ui.divV([header, list], 'pp-workspace-section');
}


export class WorkspaceTab {
  private spotlight: SpotlightWidget;

  constructor(spotlight: SpotlightWidget) {
    this.spotlight = spotlight;
  }

  build(): HTMLElement {
    const root = ui.divV([], 'pp-workspace-root');

    // Greeting — shown immediately
    const user = DG.User.current();
    const firstName = user.firstName || user.friendlyName;
    const greeting = ui.div([`Welcome back, ${firstName}.`], 'pp-workspace-greeting');
    root.appendChild(greeting);

    // Loading indicator — centered in remaining space
    const loaderHost = ui.div([ui.loader()], 'pp-workspace-loader');
    root.appendChild(loaderHost);
    this.populate(root, user).finally(() => loaderHost.remove());
    return root;
  }

  private async populate(root: HTMLElement, user: DG.User): Promise<void> {
    const [groupFavorites] = await Promise.all([
      getMyGroupFavorites(),
      //this.spotlight.initRecentData(),
    ]);

    // Group favorites sections
    if (groupFavorites.length > 0) {
      groupFavorites.sort((a, b) => a.group.friendlyName.localeCompare(b.group.friendlyName));
      for (const {group, entities, isAdmin} of groupFavorites)
        if (entities.length > 0)
          root.appendChild(createSection(group.friendlyName, 'star', entities, root, isAdmin ? group : undefined));
    }

    // New user fallback — show onboarding content instead of empty state
    if (groupFavorites.length === 0 && user.joined > dayjs().subtract(5, 'day')) {
      const onboarding = await this.spotlight.getNewUserInfoColumns();
      root.appendChild(onboarding);
    }
    else if (groupFavorites.length === 0) {
      const hint = ui.divText(
        'Your workspace is empty. Ask your admin to pin items for your group, or explore on your own.',
        'pp-workspace-empty-hint',
      );
      root.appendChild(hint);
    }

    // Recent items
    if (this.spotlight.recentEntities.length > 0) {
      const maxRecent = 5;
      root.appendChild(createRecentSection(
        this.spotlight.recentEntities.slice(0, maxRecent),
        this.spotlight.recentEntityTimes.slice(0, maxRecent),
      ));
    }
  }
}
