import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MonomerCollection, IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {Card, Control, TextInput, badge, button, component, computed, div, each, iconButton, on, show, signal, span}
  from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

import {MonomerLibManager} from './lib-manager';
import {MonomerCollectionHandler, MonomerCollectionInfo} from './monomer-collection-handler';
import {CollectionDialogs, CollectionRef} from './collection-dialogs';
import {attachMonomerCollectionsAi} from './ai-functions';
import {SEM_TYPES} from '../constants';
import {_package} from '../../package';

import '../../../css/monomer-collections.css';

/** Monomer chips a card shows before it says how many more there are. */
const CHIPS_SHOWN = 24;
const VIEW_NAME = 'Monomer Collections';

/** What the Monomer Collections view offers its ribbon, its context panel, its AI functions and its tests. */
export interface MonomerCollectionsView {
  /** The collection files, `.json` included. */
  readonly names: ReadonlySignal<string[]>;
  readonly selected: ReadonlySignal<ReadonlySet<string>>;
  readonly search: TextInput;
  /** The files the search keeps, in file order. */
  readonly shown: ReadonlySignal<string[]>;
  /** Reads the collection list and every collection again; the selection is dropped. */
  refresh(): Promise<void>;
  create(): void;
  /** Opened from the context panel's Edit. */
  editCollectionPublic(collectionName: string, collection: MonomerCollection): void;
}

/** The collection files of the package as a grid of cards: search in the ribbon, a card per
 * collection with its monomers, a New Collection card last. A click selects (Ctrl adds), a right
 * click opens the collection's menu, the context panel shows the selection. */
export function monomerCollectionsView(libHelper: MonomerLibManager, monomerLib: IMonomerLib,
  user: string): Control & MonomerCollectionsView {
  return component((self) => {
    const names = signal<string[]>([]);
    /** Each file's content once read; a file that failed to read is absent. */
    const loaded = signal<ReadonlyMap<string, MonomerCollection>>(new Map());
    const selected = signal<ReadonlySet<string>>(new Set());
    const search = new TextInput({name: 'search', search: true, width: 260,
      placeholder: 'Search by name, tag or monomer'});
    const shown = computed(() => {
      const query = search.value.trim().toLowerCase();
      const files = loaded.value;
      return names.value.filter((name) => query === '' ||
        CollectionDialogs.displayName(name).toLowerCase().includes(query) || matches(files.get(name), query));
    });
    const handler = new MonomerCollectionHandler();
    let generation = 0;

    const refresh = async (): Promise<void> => {
      const current = ++generation;
      let files: string[] = [];
      try {
        files = await libHelper.listMonomerCollections();
      }
      catch (err) {
        _package.logger.error(`Error listing monomer collections: ${err instanceof Error ? err.message : String(err)}`);
        grok.shell.error('Could not list the monomer collections');
      }
      const read = await Promise.all(files.map((name) => libHelper.readMonomerCollection(name)
        .then((data): [string, MonomerCollection] => [name, data], () => null)));
      if (current !== generation || self.isDisposed)
        return;
      select(new Set());
      loaded.value = new Map(read.filter((x): x is [string, MonomerCollection] => x !== null));
      names.value = files;
    };

    const create = (): void => {
      CollectionDialogs.edit(libHelper, monomerLib, null, () => refresh());
    };

    const editCollectionPublic = (collectionName: string, collection: MonomerCollection): void => {
      CollectionDialogs.edit(libHelper, monomerLib, {name: collectionName, data: collection}, () => refresh());
    };

    const ref = (name: string): CollectionRef | null => {
      const data = loaded.peek().get(name);
      return data === undefined ? null : {name, data};
    };

    const info = (name: string): MonomerCollectionInfo | null => {
      const r = ref(name);
      return r === null ? null : {...r, displayName: CollectionDialogs.displayName(name)};
    };

    /** The selection, and the context panel showing it — or nothing once it is empty. */
    const select = (next: Set<string>): void => {
      const had = selected.peek().size > 0;
      selected.value = next;
      const infos = [...next].map((n) => info(n)).filter((x): x is MonomerCollectionInfo => x !== null);
      if (infos.length > 0 || had)
        grok.shell.setCurrentObject(infos.length === 0 ? null : infos.length === 1 ? infos[0] : infos, false);
    };

    const pick = (name: string, add: boolean): void => {
      const next = new Set(add ? selected.peek() : []);
      if (add && next.has(name))
        next.delete(name);
      else
        next.add(name);
      select(next);
    };

    const menu = (): void => {
      const popup = DG.Menu.popup();
      const infos = [...selected.peek()].map((n) => info(n)).filter((x): x is MonomerCollectionInfo => x !== null);
      if (infos.length === 1)
        MonomerCollectionHandler.buildContextMenu(popup, infos[0], handler);
      else if (infos.length > 1)
        MonomerCollectionHandler.buildMultiContextMenu(popup, infos, handler);
      popup.show();
    };

    const tooltip = (symbol: string): HTMLElement | null => {
      const types: HelmType[] = [HelmTypes.AA, HelmTypes.NUCLEOTIDE, HelmTypes.CHEM, HelmTypes.BLOB];
      for (const type of types) {
        try {
          if (monomerLib.getWebEditorMonomer(type, symbol))
            return monomerLib.getTooltip(type, symbol);
        }
        catch {
        }
      }
      return null;
    };

    /** A monomer: its structure on hover, the monomer in the context panel on click. */
    const chip = (symbol: string): HTMLElement => {
      const el = span(symbol, 'monomer-collection-chip');
      el.setAttribute('role', 'button');
      el.tabIndex = -1;
      ui.tooltip.bind(el, () => tooltip(symbol) ?? symbol);
      return on(el, 'click', () =>
        grok.shell.setCurrentObject(DG.SemanticValue.fromValueType(symbol, SEM_TYPES.MONOMER), false));
    };

    const body = (data: MonomerCollection): HTMLElement => {
      const symbols = data.monomerSymbols ?? [];
      const chips = symbols.slice(0, CHIPS_SHOWN).map((symbol) => chip(symbol));
      if (symbols.length > CHIPS_SHOWN)
        chips.push(span(`+${symbols.length - CHIPS_SHOWN} more`, 'monomer-collection-more'));
      return div([
        ...data.description ? [span(data.description, 'monomer-collection-card-description')] : [],
        ...data.tags?.length ? [div(data.tags.map((t) => badge(t)), 'monomer-collection-card-tags')] : [],
        span(CollectionDialogs.count(symbols.length, 'monomer'), 'monomer-collection-count'),
        div(chips, 'monomer-collection-chips'),
      ], 'monomer-collection-card-content');
    };

    const actions = (name: string, data: {readonly value: MonomerCollection | undefined}): HTMLElement[] => {
      const notOwner = computed(() => data.value === undefined ? null :
        CollectionDialogs.notEditable(data.value, user));
      return [
        button('Edit', () => {
          const r = ref(name);
          if (r !== null)
            CollectionDialogs.edit(libHelper, monomerLib, r, () => refresh());
        }, {disabled: notOwner}),
        button('Delete', () => {
          const r = ref(name);
          if (r !== null)
            CollectionDialogs.delete(libHelper, [r], () => refresh());
        }, {disabled: notOwner}),
      ];
    };

    const collectionCard = (name: string): HTMLElement => {
      const title = CollectionDialogs.displayName(name);
      const data = computed(() => loaded.value.get(name));
      const isSelected = signal(false);
      const card = new Card({
        title,
        subtitle: computed(() => byline(data.value)),
        clickable: true,
        selected: isSelected,
        onClick: (e) => pick(name, e.ctrlKey || e.metaKey),
        footer: actions(name, data),
      });
      card.name = title.replace(/\s+/g, '-');
      card.root.classList.add('monomer-collection-card');
      card.body.classList.add('monomer-collection-card-body');
      // the file the card stands for, as the view before u2 marked it: lifecycle tests find a collection by it
      card.root.dataset.collectionName = name;
      card.effect(() => isSelected.value = selected.value.has(name));
      card.build(() => each(card.body, computed(() => [data.value ?? null]),
        (d) => d === null ? '' : JSON.stringify(d),
        (d) => d === null ? span('Could not read this collection', 'monomer-collection-card-error') : body(d)));
      card.scope.listen(card.root, 'contextmenu', (e) => {
        e.preventDefault();
        e.stopPropagation();
        if (!selected.peek().has(name))
          pick(name, false);
        menu();
      });
      return card.root;
    };

    const addCard = (): HTMLElement => {
      const card = new Card({icon: 'plus', title: 'New Collection', clickable: true, onClick: () => create()});
      card.name = 'New-Collection';
      card.root.classList.add('monomer-collection-add-card');
      card.body.hidden = true;
      return card.root;
    };

    const grid = div([], 'monomer-collections-grid');
    each(grid, computed(() => [...shown.value, '']), (name) => name,
      (name) => name === '' ? addCard() : collectionCard(name));
    const empty = div([], 'monomer-collections-empty');
    show(empty, computed(() => shown.value.length === 0 && names.value.length > 0),
      () => span('No collection matches the search.'));
    const root = div([empty, grid], 'monomer-collections-view');
    self.listen(root, 'click', (e) => {
      if ((e.target as Element).closest('[data-u2="card"]') === null)
        select(new Set());
    });
    self.listen(root, 'keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Escape')
        select(new Set());
    });
    self.listen(root, 'contextmenu', (e) => {
      if ((e.target as Element).closest('[data-u2="card"]') !== null)
        return;
      e.preventDefault();
      DG.Menu.popup()
        .item('New Collection...', () => create())
        .item('Refresh', () => void refresh())
        .show();
    });
    root.tabIndex = 0;
    return {root, names, selected, search, shown, refresh, create, editCollectionPublic};
  }, {name: 'monomerCollections'});
}

function matches(data: MonomerCollection | undefined, query: string): boolean {
  return data !== undefined && (data.tags?.some((t) => t.toLowerCase().includes(query)) === true ||
    data.description?.toLowerCase().includes(query) === true ||
    data.monomerSymbols?.some((s) => s.toLowerCase() === query) === true);
}

function byline(data: MonomerCollection | undefined): string {
  if (data === undefined)
    return '';
  const date = data.updatedOn ? new Date(data.updatedOn) : null;
  return [data.updatedBy ? `by ${data.updatedBy}` : '',
    date !== null && !Number.isNaN(date.getTime()) ? date.toLocaleDateString() : ''].filter((x) => x).join(' · ');
}

/** Creates and returns the Monomer Collections view; with `addView`, shows it — the open one, if any. */
export async function showMonomerCollectionsView(addView = true): Promise<DG.ViewBase> {
  if (addView) {
    const existing = Array.from(grok.shell.views).find((v) => v.name === VIEW_NAME);
    if (existing) {
      grok.shell.v = existing;
      return existing;
    }
  }
  const libHelper = await MonomerLibManager.getInstance();
  await libHelper.awaitLoaded();
  let content!: Control & MonomerCollectionsView;
  const view = ViewHost.create({name: VIEW_NAME, content: () => {
    content = monomerCollectionsView(libHelper, libHelper.getMonomerLib(), DG.User.current().login);
    return {root: content,
      ribbon: [[content.search], [
        iconButton('plus', () => content.create(), {tooltip: 'New collection'}),
        iconButton('sync', () => void content.refresh(), {tooltip: 'Refresh collections'}),
      ]],
      status: computed(() => {
        const shown = content.shown.value.length;
        const all = content.names.value.length;
        const total = CollectionDialogs.count(all, 'collection');
        return shown === all ? total : `${shown} of ${total}`;
      })};
  }}).view;
  attachMonomerCollectionsAi(view, content);
  await content.refresh();
  if (addView)
    grok.shell.addView(view);
  return view;
}
