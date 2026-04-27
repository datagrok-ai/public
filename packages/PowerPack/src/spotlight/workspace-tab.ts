import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as rxjs from 'rxjs';
import dayjs from 'dayjs';
import {getMyGroupFavorites, getAdminGroups, getMyPersonalFavorites, sortGroupsByFriendlyName, pinEntityToGroup}
  from './group-favorites';
import {showManageFavoritesDialog} from './manage-favorites-dialog';
import {SpotlightWidget} from './spotlight-widget';
import {
  clearWorkspacePreview, getWorkspacePreviewHost, showWorkspacePreview,
} from './preview-host';
import type {GroupFavorites} from './group-favorites';


/** Returns the platform icon for an entity via its registered EntityMeta handler. */
export function entityIcon(entity: DG.Entity): HTMLElement {
  return DG.ObjectHandler.forEntity(entity)?.renderIcon(entity.dart) ?? ui.iconFA('file');
}

export function isApp(entity: DG.Entity): entity is DG.Func {
  return entity instanceof DG.Func && entity.options['role'] === DG.FUNC_TYPES.APP;
}

function isRunnable(entity: DG.Entity): entity is DG.Func {
  return entity instanceof DG.Func && !isApp(entity);
}

/** Opens an entity in the appropriate way. */
function openEntity(entity: DG.Entity): void {
  if (entity instanceof DG.Project)
    entity.open();
  else if (entity instanceof DG.Func)
    entity.apply();
  else if (entity instanceof DG.FileInfo) {
    const handler = DG.ObjectHandler.forEntity(entity);
    if (handler)
      handler.renderPreview(entity).then((view) => { if (view?.root) grok.shell.addView(view); }).catch(() => { grok.shell.o = entity; });
    else
      grok.shell.o = entity;
  }
  else
    grok.shell.o = entity;
}


export class WorkspaceTab {
  private spotlight: SpotlightWidget;
  private overlayEl?: HTMLElement;
  private overlayMouseUpListener?: () => void;
  private rebuildSections?: () => Promise<void>;

  private listPane!: HTMLElement;
  private editorPane!: HTMLElement;
  private statusEl?: HTMLElement;
  private selectedRow?: HTMLElement;
  private selectedEntity?: DG.Entity;
  /** Most recent parameter editor — used to ignore results from superseded selections. */
  private activeFuncCall?: DG.FuncCall;
  /** Per-selection token — invalidates an in-flight project preview if the user moves on. */
  private activePreviewToken?: symbol;
  /** Currently embedded space view (mounted in the right pane). */
  private activeSpaceView?: DG.SpaceView;
  /** Subscription on the embedded space view's selection event. */
  private activeSpaceSub?: rxjs.Subscription;

  constructor(spotlight: SpotlightWidget) {
    this.spotlight = spotlight;
  }

  build(): HTMLElement {
    try
    {
      grok.shell.info('building');
      const root = ui.divV([], 'pp-workspace-root');

      const user = DG.User.current();
      const firstName = user.firstName || user.friendlyName;
      // const greeting = ui.div([`Welcome back, ${firstName}.`], 'pp-workspace-greeting');
      // root.appendChild(greeting);

      this.listPane = ui.divV([], 'pp-workspace-list-pane');
      this.editorPane = ui.divV([], 'pp-workspace-editor-pane');
      const main = ui.divH([this.listPane, this.editorPane], 'pp-workspace-main');
      root.appendChild(main);

      this.renderEmptyEditor();

      const loaderHost = ui.div([ui.loader()], 'pp-workspace-loader');
      this.listPane.appendChild(loaderHost);
      this.populate(root, user).finally(() => loaderHost.remove());
      return root;
    }
    catch (e) {
      console.log(e);
      return ui.divText('foo');
    }
  }

  private async populate(root: HTMLElement, user: DG.User): Promise<void> {
    const [groupFavorites, adminGroups, personalFavorites] = await Promise.all([
      getMyGroupFavorites(),
      getAdminGroups(),
      getMyPersonalFavorites(),
    ]);

    const personalContainer = ui.divV([], 'pp-workspace-personal-container');
    const groupContainer = ui.divV([], 'pp-workspace-group-container');
    const recentContainer = ui.divV([], 'pp-workspace-recent-container');
    this.listPane.append(personalContainer, groupContainer, recentContainer);

    this.rebuildSections = async (): Promise<void> => {
      const [freshFavorites, freshAdminGroups, freshPersonal] = await Promise.all([
        getMyGroupFavorites(),
        getAdminGroups(),
        getMyPersonalFavorites(),
      ]);
      personalContainer.innerHTML = '';
      groupContainer.innerHTML = '';
      this.renderPersonalSection(personalContainer, freshPersonal);
      this.renderGroupSections(groupContainer, freshFavorites, freshAdminGroups);
      this.rebuildOverlayZones(freshAdminGroups);

      // If the current selection was unpinned, clear it.
      if (this.selectedEntity && !this.listPane.querySelector(
          `.pp-workspace-item[data-entity-id="${this.selectedEntity.id}"]`))
        this.clearSelection();
    };

    this.renderPersonalSection(personalContainer, personalFavorites);
    this.renderGroupSections(groupContainer, groupFavorites, adminGroups);
    this.installDropTargets(root, adminGroups);

    if (groupFavorites.length === 0 && personalFavorites.length === 0) {
      if (user.joined > dayjs().subtract(5, 'day')) {
        const onboarding = await this.spotlight.getNewUserInfoColumns();
        this.editorPane.innerHTML = '';
        this.editorPane.appendChild(onboarding);
      }
      else {
        this.renderEmptyEditor(
          'Your workspace is empty. Drag entities here to pin them, or ask your admin to pin items for your group.',
        );
      }
    }

    // Recent items — shown at the bottom of the list pane.
    if (this.spotlight.recentEntities.length > 0) {
      const maxRecent = 5;
      recentContainer.appendChild(this.renderRecentSection(
        this.spotlight.recentEntities.slice(0, maxRecent),
        this.spotlight.recentEntityTimes.slice(0, maxRecent),
      ));
    }
  }

  private renderPersonalSection(container: HTMLElement, entities: DG.Entity[]): void {
    if (entities.length === 0)
      return;
    container.appendChild(
      this.renderSection('Myself only', 'user', entities, DG.User.current().group),
    );
  }

  private renderGroupSections(
    container: HTMLElement,
    groupFavorites: GroupFavorites[],
    adminGroups: DG.Group[],
  ): void {
    const rebuild = this.rebuildSections!;
    const displayedGroupIds = new Set<string>();

    if (groupFavorites.length > 0) {
      groupFavorites.sort((a, b) => a.group.friendlyName.localeCompare(b.group.friendlyName));
      for (const {group, entities, isAdmin} of groupFavorites) {
        if (entities.length > 0) {
          displayedGroupIds.add(group.id);
          const onEdit = isAdmin ? () => {
            showManageFavoritesDialog(group, [...entities], rebuild);
          } : undefined;
          container.appendChild(
            this.renderSection(group.friendlyName, 'star', entities, isAdmin ? group : undefined, onEdit),
          );
        }
      }
    }

    // const extraGroups = adminGroups.filter((g) => !displayedGroupIds.has(g.id));
    // if (extraGroups.length > 0) {
    //   const manageLink = ui.link('Manage favorites for another group...', () => {
    //     const menu = DG.Menu.popup();
    //     for (const g of extraGroups) {
    //       menu.item(g.friendlyName || g.name, () => {
    //         showManageFavoritesDialog(g, [], rebuild);
    //       });
    //     }
    //     menu.show();
    //   });
    //   manageLink.classList.add('pp-workspace-manage-link');
    //   container.appendChild(manageLink);
    // }
  }

  /** Section header + vertical list of entity rows. */
  private renderSection(
    title: string, iconName: string, entities: DG.Entity[],
    group?: DG.Group, onEdit?: () => void,
  ): HTMLElement {
    const headerIcon = ui.iconFA(iconName);
    const headerSpan = ui.span([headerIcon, ui.span([` ${title}`])]);
    const header = ui.h3(headerSpan, 'pp-workspace-section-header');

    if (group && onEdit) {
      const editBtn = ui.iconFA('pencil', onEdit);
      editBtn.classList.add('pp-workspace-section-edit');
      ui.tooltip.bind(editBtn, `Manage favorites for ${group.friendlyName}`);
      header.appendChild(editBtn);
    }

    const rows = ui.divV(
      entities.map((e) => this.renderItemRow(e, group)),
      'pp-workspace-items',
    );
    return ui.divV([header, rows], 'pp-workspace-section');
  }

  private renderItemRow(entity: DG.Entity, group?: DG.Group): HTMLElement {
    const icon = entityIcon(entity);
    icon.classList.add('pp-workspace-item-icon');
    const label = ui.div([entity.friendlyName], 'pp-workspace-item-label');
    const row = ui.divH([icon, label], 'pp-workspace-item');
    row.dataset.entityId = entity.id;
    ui.tooltip.bind(row, entity.friendlyName);
    ui.bind(entity, row, {contextMenu: true});

    row.addEventListener('click', () => this.selectEntity(entity, row));
    row.addEventListener('dblclick', () => openEntity(entity));

    if (group) {
      const removeBtn = ui.icons.close(async () => {
        await DG.Favorites.remove(entity, group);
        await this.rebuildSections?.();
      }, `Unpin from ${group.friendlyName}`);
      removeBtn.classList.add('pp-workspace-item-remove');
      removeBtn.addEventListener('click', (e: MouseEvent) => e.stopPropagation());
      row.appendChild(removeBtn);
    }

    return row;
  }

  private renderRecentSection(
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
      icon.classList.add('pp-workspace-item-icon');
      const label = ui.div([entity.friendlyName], 'pp-workspace-item-label');
      const row = ui.divH([icon, label], 'pp-workspace-item pp-workspace-item-recent');
      row.dataset.entityId = entity.id;

      row.addEventListener('click', () => this.selectEntity(entity, row));
      row.addEventListener('dblclick', () => openEntity(entity));
      ui.bind(entity, row, {contextMenu: true});

      const time = times[i];
      if (time) {
        const timestamp = ui.time.shortTimestamp(time);
        timestamp.classList.add('pp-workspace-item-time');
        row.appendChild(timestamp);
      }
      items.push(row);
    }
    return ui.divV([header, ui.divV(items, 'pp-workspace-items')], 'pp-workspace-section');
  }

  private selectEntity(entity: DG.Entity, row: HTMLElement): void {
    if (this.selectedEntity?.id === entity.id)
      return;
    this.selectedRow?.classList.remove('pp-workspace-item-selected');
    this.selectedRow = row;
    this.selectedEntity = entity;
    row.classList.add('pp-workspace-item-selected');
    this.renderEditorFor(entity);
  }

  private clearSelection(): void {
    this.selectedRow?.classList.remove('pp-workspace-item-selected');
    this.selectedRow = undefined;
    this.selectedEntity = undefined;
    this.activeFuncCall = undefined;
    this.activePreviewToken = undefined;
    this.statusEl = undefined;
    this.detachSpaceView();
    clearWorkspacePreview();
    this.renderEmptyEditor();
  }

  private detachSpaceView(): void {
    this.activeSpaceSub?.unsubscribe();
    this.activeSpaceSub = undefined;
    this.activeSpaceView = undefined;
  }

  private renderEmptyEditor(message?: string): void {
    this.editorPane.innerHTML = '';
    this.statusEl = undefined;
    const hint = ui.divText(
      message ?? 'Select a pinned item on the left to see its controls here.',
      'pp-workspace-editor-hint',
    );
    this.editorPane.appendChild(hint);
    clearWorkspacePreview();
  }

  private renderEditorFor(entity: DG.Entity): void {
    this.editorPane.innerHTML = '';
    this.statusEl = undefined;
    this.activePreviewToken = undefined;
    this.detachSpaceView();
    clearWorkspacePreview();

    const title = ui.divH([
      entityIcon(entity),
      ui.span([entity.friendlyName]),
    ], 'pp-workspace-editor-title');
    const openBtn = ui.button('Open', () => openEntity(entity));
    openBtn.classList.add('pp-workspace-editor-open');
    const header = ui.divH([title, openBtn], 'pp-workspace-editor-header');
    this.editorPane.appendChild(header);

    if (isApp(entity))
      this.renderAppEditor(entity);
    else if (isRunnable(entity))
      this.renderRunnableEditor(entity);
    else if (entity instanceof DG.Project && entity.isSpace)
      this.renderSpaceEditor(entity);
    else if (entity instanceof DG.Project)
      this.renderProjectEditor(entity);
    else if (entity instanceof DG.FileInfo)
      this.renderFileEditor(entity);
    else
      this.renderEntityDetails(entity);
  }

  private renderSpaceEditor(space: DG.Project): void {
    const view = DG.SpaceView.forProject(space);
    view.showItemPreview = false; // we render per-item previews in the bottom pane

    const wrap = ui.div([view.root], 'pp-workspace-space-view');
    this.editorPane.appendChild(wrap);

    this.activeSpaceView = view;
    this.activeSpaceSub = view.onCurrentObjectChanged.subscribe((item) => {
      if (this.activeSpaceView === view && item)
        this.renderItemInBottom(item);
    });
  }

  private async renderItemInBottom(item: DG.Entity): Promise<void> {
    const token = (this.activePreviewToken = Symbol('space-item-preview'));
    const title = item.friendlyName || item.name || '';
    showWorkspacePreview(
      ui.div([ui.loader()], 'pp-workspace-preview-loader'),
      title,
      () => openEntity(item),
    );

    try {
      const handler = DG.ObjectHandler.forEntity(item);
      const view: DG.View | undefined = await handler?.renderPreview(item);
      if (token !== this.activePreviewToken)
        return;
      if (view?.root instanceof HTMLElement) {
        showWorkspacePreview(
          ui.div([view.root], 'pp-workspace-preview-view'),
          title,
          () => openEntity(item),
        );
      }
      else {
        showWorkspacePreview(
          ui.divText('No preview available.', 'pp-workspace-preview-text'),
          title,
          () => openEntity(item),
        );
      }
    }
    catch (e: any) {
      if (token !== this.activePreviewToken)
        return;
      console.error('Workspace space-item preview failed', e);
      showWorkspacePreview(
        ui.divText(`Preview failed: ${e?.message ?? e}`, 'pp-workspace-preview-error'),
        title,
      );
    }
  }

  private renderProjectEditor(project: DG.Project): void {
    const handler = DG.ObjectHandler.forEntity(project);
    const card = handler?.renderCard(project);
    if (card instanceof HTMLElement) {
      const wrap = ui.div([card], 'pp-workspace-project-card');
      this.editorPane.appendChild(wrap);
    }
    else {
      this.editorPane.appendChild(ui.divText(
        'No card renderer available for this project.',
        'pp-workspace-editor-hint',
      ));
    }
    this.renderProjectPreview(project);
  }

  private async renderProjectPreview(project: DG.Project): Promise<void> {
    const token = (this.activePreviewToken = Symbol('preview'));

    const loaderContent = ui.div([ui.loader()], 'pp-workspace-preview-loader');
    showWorkspacePreview(loaderContent, project.friendlyName);

    try {
      const handler = DG.ObjectHandler.forEntity(project);
      const view: DG.View | undefined = await handler?.renderPreview(project);
      if (token !== this.activePreviewToken)
        return;
      if (view?.root instanceof HTMLElement) {
        const wrap = ui.div([view.root], 'pp-workspace-preview-project');
        showWorkspacePreview(wrap, project.friendlyName, () => project.open());
      }
      else {
        showWorkspacePreview(
          ui.divText('No preview available for this project.', 'pp-workspace-preview-text'),
          project.friendlyName,
          () => project.open(),
        );
      }
    }
    catch (e: any) {
      if (token !== this.activePreviewToken)
        return;
      console.error('Workspace project preview failed', e);
      showWorkspacePreview(
        ui.divText(`Preview failed: ${e?.message ?? e}`, 'pp-workspace-preview-error'),
        project.friendlyName,
      );
    }
  }

  private renderAppEditor(func: DG.Func): void {
    const detailsSlot = ui.div(
      func.description ? [ui.divText(func.description, 'pp-workspace-app-description')] : [],
      'pp-workspace-app-details',
    );
    this.editorPane.appendChild(detailsSlot);
    this.renderAppPreview(func, detailsSlot);
  }

  private async renderAppPreview(func: DG.Func, detailsSlot: HTMLElement): Promise<void> {
    const token = (this.activePreviewToken = Symbol('app-preview'));
    showWorkspacePreview(ui.div([ui.loader()], 'pp-workspace-preview-loader'), func.friendlyName);
    try {
      const fc = func.prepare();
      fc.adHoc = true;
      await fc.call(false, undefined, {processed: true, report: false});
      if (token !== this.activePreviewToken)
        return;
      const view = fc.getOutputParamValue();
      if (view?.root instanceof HTMLElement) {
        (<HTMLElement>view.root).classList.remove('ui-box');
        const appHeader = view.root.querySelector('.ui-app-header');
        if (appHeader) {
          appHeader.remove();
          detailsSlot.innerHTML = '';
          detailsSlot.appendChild(appHeader);
        }
        showWorkspacePreview(ui.div([view.root], 'pp-workspace-preview-app'), func.friendlyName, () => func.apply());
      }
      else
        showWorkspacePreview(ui.divText('No preview available.', 'pp-workspace-preview-text'), func.friendlyName);
    }
    catch (e: any) {
      if (token !== this.activePreviewToken)
        return;
      console.error('Workspace app preview failed', e);
      showWorkspacePreview(
        ui.divText(`Preview failed: ${e?.message ?? e}`, 'pp-workspace-preview-error'),
        func.friendlyName,
      );
    }
  }

  private renderRunnableEditor(func: DG.Func): void {
    const loader = ui.div([ui.loader()], 'pp-workspace-editor-loader');
    this.editorPane.appendChild(loader);

    const fc = func.prepare();
    fc.adHoc = true; // don't save every preview run to the user's function run history
    this.activeFuncCall = fc;

    fc.getEditor(true, true).then((editor) => {
      if (this.activeFuncCall !== fc)
        return; // superseded by a newer selection
      loader.remove();
      editor.classList.add('pp-workspace-editor-params');
      this.editorPane.appendChild(editor);

      const runBtn = ui.bigButton('Run', () => this.runAndPreview(fc));
      const statusText = ui.divText('', 'pp-workspace-editor-status');
      this.statusEl = statusText;
      const actions = ui.divH([runBtn, statusText], 'pp-workspace-editor-actions');
      this.editorPane.appendChild(actions);

      // For no-input scalars (e.g. Func returning a constant), run once immediately — not for DataQuery,
      // which may be slow or parameterized on session context.
      if (func.inputs.length === 0 && !(func instanceof DG.DataQuery))
        this.runAndPreview(fc);
    }).catch((e) => {
      if (this.activeFuncCall !== fc)
        return;
      loader.remove();
      console.error('Failed to build editor', e);
      this.editorPane.appendChild(ui.divText(
        'Could not build a parameter editor for this item.',
        'pp-workspace-editor-error',
      ));
    });
  }

  private async runAndPreview(fc: DG.FuncCall): Promise<void> {
    if (this.activeFuncCall !== fc)
      return;
    const previewHost = getWorkspacePreviewHost();
    if (!previewHost) {
      grok.shell.warning('Preview area is not available.');
      return;
    }

    if (this.statusEl) this.statusEl.textContent = 'Running...';
    const loaderContent = ui.div([ui.loader()], 'pp-workspace-preview-loader');
    showWorkspacePreview(loaderContent, fc.func.friendlyName);

    try {
      await fc.call(false, undefined, {processed: true, report: false});
      if (this.activeFuncCall !== fc)
        return; // selection changed while running
      this.renderPreviewResult(fc);
      if (this.statusEl) this.statusEl.textContent = 'Done';
    }
    catch (e: any) {
      if (this.activeFuncCall !== fc)
        return;
      console.error('Workspace preview run failed', e);
      const errorEl = ui.divText(
        `Run failed: ${fc.errorMessage ?? e?.message ?? e}`,
        'pp-workspace-preview-error',
      );
      showWorkspacePreview(errorEl, fc.func.friendlyName);
      if (this.statusEl) this.statusEl.textContent = 'Failed';
    }
  }

  private renderPreviewResult(fc: DG.FuncCall): void {
    // Prefer the function's own result views (e.g. a query's TableView with its saved layout applied).
    let views: DG.ViewBase[] = [];
    try {
      views = fc.getResultViews() ?? [];
    }
    catch (e) {
      console.warn('FuncCall.getResultViews threw', e);
    }

    if (views.length === 1) {
      const v = views[0];
      const wrap = ui.div([v.root], 'pp-workspace-preview-view');
      const title = v.name ? `${fc.func.friendlyName} — ${v.name}` : fc.func.friendlyName;
      showWorkspacePreview(wrap, title, () => grok.shell.addView(v));
      return;
    }
    if (views.length > 1) {
      const tabs: {[name: string]: HTMLElement} = {};
      for (let i = 0; i < views.length; i++) {
        const v = views[i];
        tabs[v.name || `View ${i + 1}`] = v.root;
      }
      const tabControl = ui.tabControl(tabs);
      const wrap = ui.div([tabControl.root], 'pp-workspace-preview-view');
      showWorkspacePreview(wrap, fc.func.friendlyName, () => { for (const v of views) grok.shell.addView(v); });
      return;
    }

    // No views — render the raw output value.
    const result = fc.getOutputParamValue();
    if (result instanceof DG.DataFrame) {
      const grid = DG.Viewer.grid(result);
      const gridWrap = ui.div([grid.root], 'pp-workspace-preview-grid');
      showWorkspacePreview(gridWrap, `${fc.func.friendlyName} — ${result.rowCount} rows`);
      return;
    }
    if (result instanceof HTMLElement) {
      showWorkspacePreview(result, fc.func.friendlyName);
      return;
    }
    const text = result == null ? '(no result)' : String(result);
    showWorkspacePreview(ui.divText(text, 'pp-workspace-preview-text'), fc.func.friendlyName);
  }

  private renderEntityDetails(entity: DG.Entity): void {
    const hint = ui.divText(
      'Use "Open" to work with this item.',
      'pp-workspace-editor-hint',
    );
    this.editorPane.appendChild(hint);
  }

  private renderFileEditor(file: DG.FileInfo): void {
    const rows: {[key: string]: string} = {};
    rows['Type'] = file.isDirectory ? 'Folder' : 'File';
    const path = file.fullPath || file.path;
    if (path) rows['Path'] = path;
    if (file.updatedOn) rows['Modified'] = file.updatedOn.format('YYYY-MM-DD HH:mm');
    const table = ui.tableFromMap(rows);
    table.classList.add('pp-workspace-file-details');
    this.editorPane.appendChild(table);
    this.renderFilePreview(file);
  }

  private async renderFilePreview(file: DG.FileInfo): Promise<void> {
    const token = (this.activePreviewToken = Symbol('file-preview'));
    showWorkspacePreview(ui.div([ui.loader()], 'pp-workspace-preview-loader'), file.friendlyName);
    try {
      const handler = DG.ObjectHandler.forEntity(file);
      const view: DG.View | undefined = await handler?.renderPreview(file);
      if (token !== this.activePreviewToken)
        return;
      if (view?.root instanceof HTMLElement)
        showWorkspacePreview(ui.div([view.root], 'pp-workspace-preview-file'), file.friendlyName, () => openEntity(file));
      else
        showWorkspacePreview(ui.divText('No preview available.', 'pp-workspace-preview-text'), file.friendlyName);
    }
    catch (e: any) {
      if (token !== this.activePreviewToken)
        return;
      console.error('Workspace file preview failed', e);
      showWorkspacePreview(
        ui.divText(`Preview failed: ${e?.message ?? e}`, 'pp-workspace-preview-error'),
        file.friendlyName,
      );
    }
  }

  private installDropTargets(root: HTMLElement, adminGroups: DG.Group[]): void {
    const overlay = ui.div([], 'pp-workspace-drop-overlay');
    this.overlayEl = overlay;
    root.appendChild(overlay);
    this.rebuildOverlayZones(adminGroups);

    ui.makeDroppable(root, {
      acceptDrag: (args) => {
        const o = args.dragObject;
        if (o instanceof DG.Entity && SpotlightWidget.isSpotlightEntity(o))
          this.showOverlay();
        return false;
      },
      dropIndication: false,
    });
  }

  private rebuildOverlayZones(adminGroups: DG.Group[]): void {
    if (!this.overlayEl)
      return;
    this.overlayEl.innerHTML = '';

    const makeZone = (label: string, iconName: string, group: DG.Group): void => {
      const zone = ui.divV([
        ui.div([ui.iconFA(iconName)], 'pp-workspace-drop-zone-icon'),
        ui.div([label], 'pp-workspace-drop-zone-label'),
      ], 'pp-workspace-drop-zone');
      ui.makeDroppable(zone, {
        acceptDrop: (o) =>
          o instanceof DG.Entity && SpotlightWidget.isSpotlightEntity(o),
        doDrop: async (args) => {
          const entity = args.dragObject;
          if (!(entity instanceof DG.Entity))
            return;
          try {
            const existing = await grok.dapi.entities.getFavorites(group);
            if (existing.some((e) => e.id === entity.id))
              grok.shell.info(`"${entity.friendlyName}" is already pinned to ${label}`);
            else {
              await pinEntityToGroup(entity, group);
              grok.shell.info(`Pinned "${entity.friendlyName}" to ${label}`);
              await this.rebuildSections!();
            }
          }
          catch (e) {
            console.error('Failed to pin entity', e);
            console.log(e);
            grok.shell.error(`Failed to pin "${entity.friendlyName}"`);
          }
          finally {
            this.hideOverlay();
          }
        },
        dropIndication: false,
      });
      this.overlayEl!.appendChild(zone);
    };

    makeZone('Myself only', 'user', DG.User.current().group);
    for (const g of sortGroupsByFriendlyName(adminGroups))
      makeZone(g.friendlyName, 'users', g);
  }

  private showOverlay(): void {
    if (this.overlayMouseUpListener || !this.overlayEl)
      return;
    this.overlayEl.classList.add('pp-drop-active');
    this.overlayMouseUpListener = () => this.hideOverlay();
    document.addEventListener('mouseup', this.overlayMouseUpListener, true);
  }

  private hideOverlay(): void {
    if (!this.overlayMouseUpListener || !this.overlayEl)
      return;
    this.overlayEl.classList.remove('pp-drop-active');
    document.removeEventListener('mouseup', this.overlayMouseUpListener, true);
    this.overlayMouseUpListener = undefined;
  }
}
