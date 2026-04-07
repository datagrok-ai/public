/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MonomerCollection, IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HelmType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {MonomerSelectionWidget} from '@datagrok-libraries/bio/src/utils/monomer-selection-dialog';

import {MonomerLibManager} from './lib-manager';
import {MonomerCollectionHandler, MonomerCollectionInfo} from './monomer-collection-handler';
import {SEM_TYPES} from '../constants';
import {_package} from '../../package';

//@ts-ignore
import '../../../css/monomer-collections.css';

const PAGE_SIZE = 19;

export class MonomerCollectionsView {
  private view: DG.ViewBase | null = null;
  private libHelper: MonomerLibManager | null = null;
  private monomerLib: IMonomerLib | null = null;
  private collectionNames: string[] = [];
  private filteredNames: string[] = [];
  private displayedCount = 0;
  private cardsContainer: HTMLDivElement | null = null;
  private loadMoreContainer: HTMLDivElement | null = null;
  private searchInput: DG.InputBase<string> | null = null;
  private currentUser: string = '';
  /** Cached collection data for search filtering */
  private collectionsCache: Map<string, MonomerCollection> = new Map();
  /** Currently selected collection names */
  private selectedCollections: Set<string> = new Set();
  /** Map of collection name -> card element for quick access */
  private cardElements: Map<string, HTMLDivElement> = new Map();
  /** The handler instance for context menu and actions */
  private handler: MonomerCollectionHandler = new MonomerCollectionHandler();

  static readonly VIEW_NAME = 'Monomer Collections';

  async init(): Promise<DG.ViewBase> {
    this.libHelper = await MonomerLibManager.getInstance();
    await this.libHelper.awaitLoaded();
    this.monomerLib = this.libHelper.getMonomerLib();
    this.currentUser = DG.User.current().login;

    this.view = DG.View.create();
    this.view.name = MonomerCollectionsView.VIEW_NAME;
    // Store reference for the handler to find
    (this.view as any)._mcView = this;

    // Ribbon
    const addBtn = ui.icons.add(() => this.showAddCollectionDialog(), 'Create new monomer collection');
    const refreshBtn = ui.iconFA('sync', () => this.refresh(), 'Refresh collections');
    this.view.setRibbonPanels([[addBtn, refreshBtn]]);

    // Background context menu (not on cards)
    this.view.root.addEventListener('contextmenu', (e: MouseEvent) => {
      const target = e.target as HTMLElement;
      if (target.closest('.monomer-collection-card') || target.closest('.monomer-collection-add-card'))
        return;
      e.preventDefault();
      const menu = DG.Menu.popup();
      menu.item('New Collection...', () => this.showAddCollectionDialog());
      menu.item('Refresh', () => this.refresh());
      menu.show();
    });

    // Click on empty space deselects all
    this.view.root.addEventListener('click', (e: MouseEvent) => {
      const target = e.target as HTMLElement;
      if (!target.closest('.monomer-collection-card') && !target.closest('.monomer-collection-add-card'))
        this.clearSelection();
    });

    // Escape key clears selection
    this.view.root.tabIndex = 0;
    this.view.root.addEventListener('keydown', (e: KeyboardEvent) => {
      if (e.key === 'Escape')
        this.clearSelection();
    });

    await this.buildContent();
    return this.view;
  }

  private async buildContent(): Promise<void> {
    if (!this.view) return;
    this.view.root.innerHTML = '';

    const root = ui.div([], {classes: 'monomer-collections-view'});

    // Search bar
    this.searchInput = ui.input.string('Monomer Collections', {value: '', placeholder: 'Search by name or tag...'});
    this.searchInput.root.classList.add('monomer-collections-search');
    const searchEl = this.searchInput.input as HTMLInputElement;
    searchEl.style.width = '100%';
    searchEl.style.marginBottom = '12px';
    let searchTimeout: ReturnType<typeof setTimeout> | null = null;
    searchEl.addEventListener('input', () => {
      if (searchTimeout) clearTimeout(searchTimeout);
      searchTimeout = setTimeout(() => this.applyFilter(), 250);
    });

    this.cardsContainer = ui.div([], {classes: 'monomer-collections-grid'}) as HTMLDivElement;
    this.loadMoreContainer = ui.div([], {classes: 'monomer-collections-load-more'}) as HTMLDivElement;

    root.appendChild(this.searchInput.root);
    root.appendChild(this.cardsContainer);
    root.appendChild(this.loadMoreContainer);
    this.view.root.appendChild(root);

    await this.loadCollections();
  }

  private async loadCollections(): Promise<void> {
    try {
      this.collectionNames = await this.libHelper!.listMonomerCollections();
    } catch (err: any) {
      _package.logger.error(`Error listing monomer collections: ${err instanceof Error ? err.message : err.toString()}`);
      grok.shell.error('Error loading monomer collections');
      this.collectionNames = [];
    }

    // Preload all collection metadata in the background for search filtering
    this.collectionsCache.clear();
    for (const name of this.collectionNames) {
      this.libHelper!.readMonomerCollection(name).then((c) => {
        this.collectionsCache.set(name, c);
      }).catch(() => { /* ignore load errors for cache */ });
    }

    this.applyFilter();
  }

  private applyFilter(): void {
    const query = (this.searchInput?.value ?? '').trim().toLowerCase();
    this.filteredNames = !query ? [...this.collectionNames] :
      this.collectionNames.filter((name) => {
        const displayName = name.replace(/\.json$/i, '').toLowerCase();
        if (displayName.includes(query)) return true;
        const cached = this.collectionsCache.get(name);
        if (cached?.tags?.some((t) => t.toLowerCase().includes(query)))
          return true;
        return false;
      });

    this.displayedCount = 0;
    this.cardsContainer!.innerHTML = '';
    this.cardElements.clear();
    this.showNextPage();
  }

  private showNextPage(): void {
    const end = Math.min(this.displayedCount + PAGE_SIZE, this.filteredNames.length);
    for (let i = this.displayedCount; i < end; i++) {
      const card = this.createCollectionCard(this.filteredNames[i]);
      this.cardsContainer!.appendChild(card);
      this.cardElements.set(this.filteredNames[i], card);
    }

    this.displayedCount = end;

    // Add the "add new" card at the end (remove/re-add so it stays last)
    const existingAddCard = this.cardsContainer!.querySelector('.monomer-collection-add-card');
    if (existingAddCard) existingAddCard.remove();
    this.cardsContainer!.appendChild(this.createAddCard());

    // Show empty state if no collections
    const existingEmpty = this.cardsContainer!.querySelector('.monomer-collection-empty-state');
    if (existingEmpty) existingEmpty.remove();
    if (this.filteredNames.length === 0) {
      const msg = this.collectionNames.length === 0 ?
        'No monomer collections found. Click "+" to create one.' :
        'No collections match your search.';
      const emptyState = ui.div([ui.divText(msg)], {classes: 'monomer-collection-empty-state'});
      this.cardsContainer!.prepend(emptyState);
    }

    // Load more button
    this.loadMoreContainer!.innerHTML = '';
    if (this.displayedCount < this.filteredNames.length) {
      const remaining = this.filteredNames.length - this.displayedCount;
      const loadMoreBtn = ui.button(`Load more (${remaining} remaining)`, () => this.showNextPage());
      this.loadMoreContainer!.appendChild(loadMoreBtn);
    }
  }

  private createCollectionCard(collectionName: string): HTMLDivElement {
    const displayName = collectionName.replace(/\.json$/i, '');

    // Header with name shown immediately
    const titleEl = ui.div([displayName], {classes: 'monomer-collection-card-title'});
    ui.tooltip.bind(titleEl, displayName);
    const headerEl = ui.div([titleEl], {classes: 'monomer-collection-card-header'});

    // Body: loaded asynchronously
    const bodyEl = ui.div([], {classes: 'monomer-collection-card-body'});

    // Actions placeholder
    const actionsEl = ui.div([], {classes: 'monomer-collection-card-actions'});

    const card = ui.div([headerEl, bodyEl, actionsEl], {classes: 'monomer-collection-card'}) as HTMLDivElement;
    card.dataset.collectionName = collectionName;

    // Apply selected state if already selected
    if (this.selectedCollections.has(collectionName))
      card.classList.add('monomer-collection-card-selected');

    // Click handler: select card and set as current object
    card.addEventListener('click', (e: MouseEvent) => {
      // Don't handle clicks on buttons or monomer tags
      const target = e.target as HTMLElement;
      if (target.closest('.monomer-collection-card-actions') || target.closest('.monomer-collection-tag'))
        return;

      e.stopPropagation();
      this.handleCardClick(collectionName, e.ctrlKey || e.metaKey);
    });

    // Right-click context menu
    card.addEventListener('contextmenu', (e: MouseEvent) => {
      e.preventDefault();
      e.stopPropagation();

      // If right-clicking a non-selected card, select it first
      if (!this.selectedCollections.has(collectionName))
        this.handleCardClick(collectionName, false);

      this.showCardContextMenu(e);
    });

    // Load content in background
    const loadContent = async (): Promise<HTMLElement> => {
      const collection = await this.libHelper!.readMonomerCollection(collectionName);
      // Cache for search
      this.collectionsCache.set(collectionName, collection);

      // Description
      if (collection.description) {
        const descEl = ui.div([collection.description], {classes: 'monomer-collection-card-description'});
        ui.tooltip.bind(descEl, collection.description);
        headerEl.appendChild(descEl);
      }

      // Collection tags (not monomer symbols - these are user-defined labels)
      if (collection.tags && collection.tags.length > 0) {
        const collTagsContainer = ui.div([], {classes: 'monomer-collection-card-tags'});
        for (const t of collection.tags) {
          const tagEl = ui.div([t], {classes: 'monomer-collection-card-tag'});
          collTagsContainer.appendChild(tagEl);
        }
        headerEl.appendChild(collTagsContainer);
      }

      // Meta info
      const metaParts: string[] = [];
      if (collection.updatedBy) metaParts.push(`by ${collection.updatedBy}`);
      if (collection.updatedOn)
        try { metaParts.push(new Date(collection.updatedOn).toLocaleDateString()); } catch { /* ignore */ }
      if (metaParts.length > 0)
        headerEl.appendChild(ui.div([metaParts.join(' | ')], {classes: 'monomer-collection-card-meta'}));

      // Monomer tags
      const monomerTagsContainer = ui.div([], {classes: 'monomer-collection-tags'});
      const symbols = collection.monomerSymbols ?? [];
      for (const symbol of symbols) {
        const tag = ui.div([symbol], {classes: 'monomer-collection-tag'});
        tag.addEventListener('mouseenter', () => {
          if (this.monomerLib) {
            const tooltipEl = this.getMonomerTooltipSafe(symbol);
            if (tooltipEl) {
              const rect = tag.getBoundingClientRect();
              ui.tooltip.show(tooltipEl, rect.left, rect.bottom + 4);
            }
          }
        });
        tag.addEventListener('mouseleave', () => ui.tooltip.hide());
        // Clicking a monomer tag sets its SemanticValue as the current object
        tag.addEventListener('click', (e: MouseEvent) => {
          e.stopPropagation();
          grok.shell.o = DG.SemanticValue.fromValueType(symbol, SEM_TYPES.MONOMER);
        });
        monomerTagsContainer.appendChild(tag);
      }

      const countLabel = ui.divText(`${symbols.length} monomer(s)`, {style: {fontSize: '11px', color: 'var(--grey-4)', marginBottom: '6px'}});
      const contentDiv = ui.divV([countLabel, monomerTagsContainer]);

      // Actions: show edit/delete with disabled state for non-owners
      const isOwner = collection.updatedBy === this.currentUser;
      const editBtn = ui.button('Edit', () => { if (isOwner) this.showEditCollectionDialog(collectionName, collection); });
      const deleteBtn = ui.button('Delete', () => { if (isOwner) this.confirmDeleteCollection(collectionName); });
      if (!isOwner) {
        for (const btn of [editBtn, deleteBtn]) {
          btn.style.opacity = '0.4';
          btn.style.cursor = 'default';
        }
        ui.tooltip.bind(editBtn, 'Only the author can edit');
        ui.tooltip.bind(deleteBtn, 'Only the author can delete');
      }
      actionsEl.appendChild(editBtn);
      actionsEl.appendChild(deleteBtn);

      return contentDiv;
    };

    bodyEl.appendChild(ui.wait(() => loadContent()));
    return card;
  }

  /** Handle card click with optional Ctrl/Cmd for multi-select. */
  private handleCardClick(collectionName: string, isMulti: boolean): void {
    if (isMulti) {
      // Toggle selection
      if (this.selectedCollections.has(collectionName))
        this.selectedCollections.delete(collectionName);
      else
        this.selectedCollections.add(collectionName);
    } else {
      // Single select: clear others
      this.selectedCollections.clear();
      this.selectedCollections.add(collectionName);
    }

    this.updateCardSelectionStyles();
    this.setCurrentObject();
  }

  /** Clear all selections. */
  private clearSelection(): void {
    this.selectedCollections.clear();
    this.updateCardSelectionStyles();
  }

  /** Update card CSS classes to reflect current selection. */
  private updateCardSelectionStyles(): void {
    for (const [name, card] of this.cardElements) {
      if (this.selectedCollections.has(name))
        card.classList.add('monomer-collection-card-selected');
      else
        card.classList.remove('monomer-collection-card-selected');
    }
  }

  /** Set the current object in grok.shell.o based on selection. */
  private async setCurrentObject(): Promise<void> {
    const selected = [...this.selectedCollections];
    if (selected.length === 0) return;

    // Build MonomerCollectionInfo items for the selected collections
    const infos: MonomerCollectionInfo[] = [];
    for (const name of selected) {
      let data = this.collectionsCache.get(name);
      if (!data) {
        try { data = await this.libHelper!.readMonomerCollection(name); } catch { continue; }
      }
      infos.push({
        name,
        displayName: name.replace(/\.json$/i, ''),
        data,
      });
    }

    if (infos.length === 1)
      grok.shell.o = infos[0];
    else if (infos.length > 1)
      grok.shell.o = infos;
  }

  /** Show context menu for selected cards. */
  private showCardContextMenu(e: MouseEvent): void {
    const menu = DG.Menu.popup();
    const selected = [...this.selectedCollections];

    if (selected.length <= 1 && selected.length > 0) {
      const name = selected[0];
      const data = this.collectionsCache.get(name);
      if (data) {
        const info: MonomerCollectionInfo = {name, displayName: name.replace(/\.json$/i, ''), data};
        MonomerCollectionHandler.buildContextMenu(menu, info, this.handler);
      }
    } else if (selected.length > 1) {
      const items: MonomerCollectionInfo[] = selected
        .map((name) => {
          const data = this.collectionsCache.get(name);
          return data ? {name, displayName: name.replace(/\.json$/i, ''), data} : null;
        })
        .filter((x): x is MonomerCollectionInfo => x !== null);
      if (items.length > 0)
        MonomerCollectionHandler.buildMultiContextMenu(menu, items, this.handler);
    }

    menu.show();
  }

  private getMonomerTooltipSafe(symbol: string): HTMLElement | null {
    if (!this.monomerLib) return null;
    const helmTypes: HelmType[] = [HelmTypes.AA, HelmTypes.NUCLEOTIDE, HelmTypes.CHEM, HelmTypes.BLOB];
    for (const ht of helmTypes) {
      try {
        const wem = this.monomerLib.getWebEditorMonomer(ht, symbol);
        if (wem) return this.monomerLib.getTooltip(ht, symbol);
      } catch (_) { /* skip */ }
    }
    return null;
  }

  private createAddCard(): HTMLDivElement {
    const icon = ui.iconFA('plus', () => {}, '');
    icon.style.fontSize = '32px';
    icon.style.marginBottom = '8px';
    const label = ui.span(['New Collection']);

    const content = ui.div([icon, label], {classes: 'monomer-collection-add-card-content'});
    const card = ui.div([content], {classes: 'monomer-collection-add-card'}) as HTMLDivElement;
    card.addEventListener('click', () => this.showAddCollectionDialog());
    return card;
  }

  private showAddCollectionDialog(): void {
    const nameInput = ui.input.string('Name', {nullable: false, placeholder: 'Enter collection name'});
    const descInput = ui.input.string('Description', {nullable: true, placeholder: 'Optional description'});
    const tagsInput = ui.input.string('Tags', {nullable: true, placeholder: 'Comma-separated tags'});

    const polymerTypes: PolymerType[] = ['PEPTIDE', 'RNA', 'CHEM'];
    const polymerTypeInput = ui.input.choice('Polymer Type', {items: polymerTypes, value: 'PEPTIDE' as PolymerType, nullable: false});

    const monomerWidget = new MonomerSelectionWidget(
      this.monomerLib!, polymerTypeInput.value as PolymerType, [], this.libHelper!,
    );

    // Update widget when polymer type changes
    // eslint-disable-next-line rxjs/no-ignored-subscription
    polymerTypeInput.onChanged.subscribe(() => {
      monomerWidget.setPolymerType(polymerTypeInput.value as PolymerType);
    });

    const dlg = ui.dialog({title: 'New Monomer Collection'})
      .add(ui.divV([
        nameInput,
        descInput,
        tagsInput,
        polymerTypeInput,
        ui.element('hr'),
        monomerWidget.root,
      ]))
      .onOK(async () => {
        const name = nameInput.value?.trim();
        if (!name) {
          grok.shell.warning('Collection name is required');
          return;
        }
        const selectedMonomers = monomerWidget.getSelectedMonomers();
        if (selectedMonomers.length === 0) {
          grok.shell.warning('Please select at least one monomer');
          return;
        }
        const tags = parseTags(tagsInput.value);
        try {
          await this.libHelper!.addOrUpdateMonomerCollection(name, selectedMonomers, descInput.value ?? undefined, tags.length > 0 ? tags : undefined);
          grok.shell.info(`Collection "${name}" created successfully`);
          await this.refresh();
        } catch (err: any) {
          grok.shell.error('Error creating collection');
          _package.logger.error(`Error creating collection: ${err instanceof Error ? err.message : err.toString()}`);
        }
      })
      .show({resizable: true});
    dlg.root.style.width = '550px';
    monomerWidget.focus();
  }

  /** Public method so the handler can trigger edit from the context panel. */
  editCollectionPublic(collectionName: string, collection: MonomerCollection): void {
    this.showEditCollectionDialog(collectionName, collection);
  }

  private showEditCollectionDialog(collectionName: string, collection: MonomerCollection): void {
    const displayName = collectionName.replace(/\.json$/i, '');
    const nameInput = ui.input.string('Name', {value: displayName, nullable: false});
    const descInput = ui.input.string('Description', {value: collection.description ?? '', nullable: true});
    const tagsInput = ui.input.string('Tags', {value: (collection.tags ?? []).join(', '), nullable: true, placeholder: 'Comma-separated tags'});

    const polymerTypes: PolymerType[] = ['PEPTIDE', 'RNA', 'CHEM'];
    const polymerTypeInput = ui.input.choice('Polymer Type', {items: polymerTypes, value: 'PEPTIDE' as PolymerType, nullable: false});

    const monomerWidget = new MonomerSelectionWidget(
      this.monomerLib!, polymerTypeInput.value as PolymerType, [...(collection.monomerSymbols ?? [])], this.libHelper!,
    );

    // eslint-disable-next-line rxjs/no-ignored-subscription
    polymerTypeInput.onChanged.subscribe(() => {
      monomerWidget.setPolymerType(polymerTypeInput.value as PolymerType);
    });

    const dlg = ui.dialog({title: `Edit Collection: ${displayName}`})
      .add(ui.divV([
        nameInput,
        descInput,
        tagsInput,
        polymerTypeInput,
        ui.element('hr'),
        monomerWidget.root,
      ]))
      .onOK(async () => {
        const name = nameInput.value?.trim();
        if (!name) {
          grok.shell.warning('Collection name is required');
          return;
        }
        const selectedMonomers = monomerWidget.getSelectedMonomers();
        if (selectedMonomers.length === 0) {
          grok.shell.warning('Please select at least one monomer');
          return;
        }
        const tags = parseTags(tagsInput.value);
        try {
          if (name !== displayName)
            await this.libHelper!.deleteMonomerCollection(collectionName);

          await this.libHelper!.addOrUpdateMonomerCollection(name, selectedMonomers, descInput.value ?? undefined, tags.length > 0 ? tags : undefined);
          grok.shell.info(`Collection "${name}" saved successfully`);
          await this.refresh();
        } catch (err: any) {
          grok.shell.error('Error saving collection');
          _package.logger.error(`Error saving collection: ${err instanceof Error ? err.message : err.toString()}`);
        }
      })
      .show({resizable: true});
    dlg.root.style.width = '550px';
    monomerWidget.focus();
  }

  private confirmDeleteCollection(collectionName: string): void {
    const displayName = collectionName.replace(/\.json$/i, '');
    ui.dialog({title: 'Delete Collection'})
      .add(ui.divText(`Are you sure you want to delete the collection "${displayName}"?`))
      .onOK(async () => {
        try {
          await this.libHelper!.deleteMonomerCollection(collectionName);
          grok.shell.info(`Collection "${displayName}" deleted`);
          await this.refresh();
        } catch (err: any) {
          grok.shell.error('Error deleting collection');
          _package.logger.error(`Error deleting collection: ${err instanceof Error ? err.message : err.toString()}`);
        }
      })
      .show();
  }

  async refresh(): Promise<void> {
    this.selectedCollections.clear();
    await this.loadCollections();
  }
}

function parseTags(input: string | null): string[] {
  if (!input) return [];
  return input.split(',').map((t) => t.trim()).filter((t) => t.length > 0);
}

/** Creates and returns the Monomer Collections management view */
export async function showMonomerCollectionsView(addView = true): Promise<DG.ViewBase> {
  if (addView) {
    const existingView = Array.from(grok.shell.views).find(
      (v) => v.name === MonomerCollectionsView.VIEW_NAME,
    );
    if (existingView) {
      grok.shell.v = existingView;
      return existingView;
    }
  }

  const manager = new MonomerCollectionsView();
  const view = await manager.init();
  if (addView)
    grok.shell.addView(view);
  return view;
}
