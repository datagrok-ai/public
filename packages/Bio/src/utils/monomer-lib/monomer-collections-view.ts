import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MonomerCollection, IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {Card, Control, TextInput, badge, button, computed, div, each, iconButton, on, show, signal, span}
  from '@datagrok-libraries/u2';
import {appView} from '@datagrok-libraries/u2/src/dg/index.js';

import {MonomerLibManager} from './lib-manager';
import {MonomerCollectionHandler, MonomerCollectionInfo} from './monomer-collection-handler';
import {CollectionDialogs, CollectionRef} from './collection-dialogs';
import {attachMonomerCollectionsAi} from './ai-functions';
import {SEM_TYPES} from '../constants';
import {_package} from '../../package';

import '../../../css/monomer-collections.css';

/** Monomer chips a card shows before it says how many more there are. */
const CHIPS_SHOWN = 24;

/** The collection files of the package as a grid of cards: search in the ribbon, a card per
 * collection with its monomers, a New Collection card last. A click selects (Ctrl adds), a right
 * click opens the collection's menu, the context panel shows the selection. */
export class MonomerCollectionsView extends Control {
  static readonly VIEW_NAME = 'Monomer Collections';

  /** The collection files, `.json` included. */
  readonly names = signal<string[]>([]);
  /** Each file's content once read; a file that failed to read is absent. */
  readonly loaded = signal<ReadonlyMap<string, MonomerCollection>>(new Map());
  readonly selected = signal<ReadonlySet<string>>(new Set());
  readonly search = new TextInput({name: 'search', search: true, width: 260,
    placeholder: 'Search by name, tag or monomer'});
  /** The files the search keeps, in file order. */
  readonly shown = computed(() => {
    const query = this.search.value.value.trim().toLowerCase();
    const loaded = this.loaded.value;
    return this.names.value.filter((name) => query === '' ||
      CollectionDialogs.displayName(name).toLowerCase().includes(query) ||
      MonomerCollectionsView._matches(loaded.get(name), query));
  });

  private readonly _handler = new MonomerCollectionHandler();
  private _generation = 0;

  /** The Monomer Collections view open in the shell, if any. */
  static current(): MonomerCollectionsView | null {
    const view = Array.from(grok.shell.views).find((v) => v.name === MonomerCollectionsView.VIEW_NAME);
    const root = view?.root.querySelector<HTMLElement>('.monomer-collections-view');
    const content = root ? Control.forElement(root) : undefined;
    return content instanceof MonomerCollectionsView ? content : null;
  }

  constructor(private readonly libHelper: MonomerLibManager, private readonly monomerLib: IMonomerLib,
    private readonly user: string) {
    super();
    this.root.classList.add('monomer-collections-view');
    this.name = 'monomerCollections';
    // built with the fields, before any owner: it lives in the ribbon but goes with the view
    this.own(() => this.search.dispose());
    this.runInScope(() => {
      const grid = div([], 'monomer-collections-grid');
      each(grid, computed(() => [...this.shown.value, '']), (name) => name,
        (name) => name === '' ? this._addCard() : this._card(name));
      const empty = div([], 'monomer-collections-empty');
      show(empty, computed(() => this.shown.value.length === 0 && this.names.value.length > 0),
        () => span('No collection matches the search.'));
      this.root.append(empty, grid);
    });
    this.scope.listen(this.root, 'click', (e) => {
      if ((e.target as Element).closest('[data-u2="card"]') === null)
        this._select(new Set());
    });
    this.scope.listen(this.root, 'keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Escape')
        this._select(new Set());
    });
    this.scope.listen(this.root, 'contextmenu', (e) => {
      if ((e.target as Element).closest('[data-u2="card"]') !== null)
        return;
      e.preventDefault();
      DG.Menu.popup()
        .item('New Collection...', () => this.create())
        .item('Refresh', () => void this.refresh())
        .show();
    });
    this.root.tabIndex = 0;
  }

  /** Reads the collection list and every collection again; the selection is dropped. */
  async refresh(): Promise<void> {
    const generation = ++this._generation;
    let names: string[] = [];
    try {
      names = await this.libHelper.listMonomerCollections();
    } catch (err) {
      _package.logger.error(`Error listing monomer collections: ${err instanceof Error ? err.message : String(err)}`);
      grok.shell.error('Could not list the monomer collections');
    }
    const read = await Promise.all(names.map((name) => this.libHelper.readMonomerCollection(name)
      .then((data): [string, MonomerCollection] => [name, data], () => null)));
    if (generation !== this._generation || this.scope.isDisposed)
      return;
    this._select(new Set());
    this.loaded.value = new Map(read.filter((x): x is [string, MonomerCollection] => x !== null));
    this.names.value = names;
  }

  create(): void {
    CollectionDialogs.edit(this.libHelper, this.monomerLib, null, () => this.refresh());
  }

  /** Opened from the context panel's Edit. */
  editCollectionPublic(collectionName: string, collection: MonomerCollection): void {
    CollectionDialogs.edit(this.libHelper, this.monomerLib, {name: collectionName, data: collection},
      () => this.refresh());
  }

  private _ref(name: string): CollectionRef | null {
    const data = this.loaded.peek().get(name);
    return data === undefined ? null : {name, data};
  }

  private _info(name: string): MonomerCollectionInfo | null {
    const ref = this._ref(name);
    return ref === null ? null : {...ref, displayName: CollectionDialogs.displayName(name)};
  }

  private _card(name: string): HTMLElement {
    const title = CollectionDialogs.displayName(name);
    const data = computed(() => this.loaded.value.get(name));
    const isSelected = signal(false);
    const card = new Card({
      title,
      subtitle: computed(() => MonomerCollectionsView._byline(data.value)),
      clickable: true,
      selected: isSelected,
      onClick: (e) => this._pick(name, e.ctrlKey || e.metaKey),
      footer: this._actions(name, data),
    });
    card.name = title.replace(/\s+/g, '-');
    card.root.classList.add('monomer-collection-card');
    // the file the card stands for, as the view before u2 marked it: lifecycle tests find a collection by it
    card.root.dataset.collectionName = name;
    card.effect(() => isSelected.value = this.selected.value.has(name));
    card.runInScope(() => each(card.body, computed(() => [data.value ?? null]),
      (d) => d === null ? '' : JSON.stringify(d),
      (d) => d === null ? span('Could not read this collection', 'monomer-collection-card-error') : this._body(d)));
    card.scope.listen(card.root, 'contextmenu', (e) => {
      e.preventDefault();
      e.stopPropagation();
      if (!this.selected.peek().has(name))
        this._pick(name, false);
      this._menu();
    });
    return card.root;
  }

  private _actions(name: string, data: {readonly value: MonomerCollection | undefined}): HTMLElement[] {
    const notOwner = computed(() => data.value === undefined ? null :
      CollectionDialogs.notEditable(data.value, this.user));
    return [
      button('Edit', () => {
        const ref = this._ref(name);
        if (ref !== null)
          CollectionDialogs.edit(this.libHelper, this.monomerLib, ref, () => this.refresh());
      }, {disabled: notOwner}),
      button('Delete', () => {
        const ref = this._ref(name);
        if (ref !== null)
          CollectionDialogs.delete(this.libHelper, [ref], () => this.refresh());
      }, {disabled: notOwner}),
    ];
  }

  private _body(data: MonomerCollection): HTMLElement {
    const symbols = data.monomerSymbols ?? [];
    const chips = symbols.slice(0, CHIPS_SHOWN).map((symbol) => this._chip(symbol));
    if (symbols.length > CHIPS_SHOWN)
      chips.push(span(`+${symbols.length - CHIPS_SHOWN} more`, 'monomer-collection-more'));
    return div([
      ...data.description ? [span(data.description, 'monomer-collection-card-description')] : [],
      ...data.tags?.length ? [div(data.tags.map((t) => badge(t)), 'monomer-collection-card-tags')] : [],
      span(CollectionDialogs.count(symbols.length, 'monomer'), 'monomer-collection-count'),
      div(chips, 'monomer-collection-chips'),
    ], 'monomer-collection-card-content');
  }

  /** A monomer: its structure on hover, the monomer in the context panel on click. */
  private _chip(symbol: string): HTMLElement {
    const chip = span(symbol, 'monomer-collection-chip');
    chip.setAttribute('role', 'button');
    chip.tabIndex = -1;
    ui.tooltip.bind(chip, () => this._tooltip(symbol) ?? symbol);
    return on(chip, 'click', () =>
      grok.shell.setCurrentObject(DG.SemanticValue.fromValueType(symbol, SEM_TYPES.MONOMER), false));
  }

  private _addCard(): HTMLElement {
    const card = new Card({icon: 'plus', title: 'New Collection', clickable: true, onClick: () => this.create()});
    card.name = 'New-Collection';
    card.root.classList.add('monomer-collection-add-card');
    return card.root;
  }

  private _pick(name: string, add: boolean): void {
    const next = new Set(add ? this.selected.peek() : []);
    if (add && next.has(name))
      next.delete(name);
    else
      next.add(name);
    this._select(next);
  }

  /** The selection, and the context panel showing it — or nothing once it is empty. */
  private _select(names: Set<string>): void {
    const had = this.selected.peek().size > 0;
    this.selected.value = names;
    const infos = [...names].map((n) => this._info(n)).filter((x): x is MonomerCollectionInfo => x !== null);
    if (infos.length > 0 || had)
      grok.shell.setCurrentObject(infos.length === 0 ? null : infos.length === 1 ? infos[0] : infos, false);
  }

  private _menu(): void {
    const menu = DG.Menu.popup();
    const infos = [...this.selected.peek()].map((n) => this._info(n))
      .filter((x): x is MonomerCollectionInfo => x !== null);
    if (infos.length === 1)
      MonomerCollectionHandler.buildContextMenu(menu, infos[0], this._handler);
    else if (infos.length > 1)
      MonomerCollectionHandler.buildMultiContextMenu(menu, infos, this._handler);
    menu.show();
  }

  private _tooltip(symbol: string): HTMLElement | null {
    const types: HelmType[] = [HelmTypes.AA, HelmTypes.NUCLEOTIDE, HelmTypes.CHEM, HelmTypes.BLOB];
    for (const type of types) {
      try {
        if (this.monomerLib.getWebEditorMonomer(type, symbol))
          return this.monomerLib.getTooltip(type, symbol);
      } catch {
      }
    }
    return null;
  }

  private static _matches(data: MonomerCollection | undefined, query: string): boolean {
    return data !== undefined && (data.tags?.some((t) => t.toLowerCase().includes(query)) === true ||
      data.description?.toLowerCase().includes(query) === true ||
      data.monomerSymbols?.some((s) => s.toLowerCase() === query) === true);
  }

  private static _byline(data: MonomerCollection | undefined): string {
    if (data === undefined)
      return '';
    const date = data.updatedOn ? new Date(data.updatedOn) : null;
    return [data.updatedBy ? `by ${data.updatedBy}` : '',
      date !== null && !Number.isNaN(date.getTime()) ? date.toLocaleDateString() : ''].filter((x) => x).join(' · ');
  }
}

/** Creates and returns the Monomer Collections view; with `addView`, shows it — the open one, if any. */
export async function showMonomerCollectionsView(addView = true): Promise<DG.ViewBase> {
  if (addView) {
    const existing = Array.from(grok.shell.views).find((v) => v.name === MonomerCollectionsView.VIEW_NAME);
    if (existing) {
      grok.shell.v = existing;
      return existing;
    }
  }
  const libHelper = await MonomerLibManager.getInstance();
  await libHelper.awaitLoaded();
  const content = new MonomerCollectionsView(libHelper, libHelper.getMonomerLib(), DG.User.current().login);
  const view = appView({
    name: MonomerCollectionsView.VIEW_NAME,
    content,
    ribbon: [[content.search], [
      iconButton('plus', () => content.create(), {tooltip: 'New collection'}),
      iconButton('sync', () => void content.refresh(), {tooltip: 'Refresh collections'}),
    ]],
    status: computed(() => {
      const shown = content.shown.value.length;
      const all = content.names.value.length;
      const total = CollectionDialogs.count(all, 'collection');
      return shown === all ? total : `${shown} of ${total}`;
    }),
  });
  attachMonomerCollectionsAi(view, content);
  await content.refresh();
  if (addView)
    grok.shell.addView(view);
  return view;
}
