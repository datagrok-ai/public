import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Control} from '../../core/component.js';
import {Input} from '../../core/input-base.js';
import {computed, signal, ReadonlySignal} from '../../core/signals.js';
import {div, divH, span} from '../../core/elements.js';
import {icon} from '../../components/display/icon.js';
import {Dialog} from '../../components/containers/dialog.js';
import {TextInput, TextArea} from '../../components/inputs/text-input.js';
import {VirtualTree} from '../../components/collections/tree.js';
import type {TreeNode} from '../../components/collections/tree.js';
import {VirtualList} from '../../components/collections/list.js';
import {Format} from '../../core/format.js';
import {entityCard, handlerRenderer} from '../entities/entity.js';
import {SpaceInput} from '../inputs/space-input.js';
import {SpaceTree, SpaceEntry, spaceHasName, rememberSpace} from './spaces.js';

export interface SaveToSpaceValues {
  name: string;
  description: string;
}

export interface SaveToSpaceOptions {
  title?: string;
  name?: string;
  /** The description box; left out without one. */
  description?: string;
  space?: DG.Project | null;
  /** Offers this space and the spaces under it alone. */
  root?: DG.Project;
  /** More inputs, under the name and the description. */
  inputs?: Input<any>[];
  /** Whether `space` holds `name` already; `subspaceExists` of the space by default. */
  exists?: (space: DG.Project, name: string) => Promise<boolean>;
  /** The save itself. */
  save: (space: DG.Project, values: SaveToSpaceValues) => unknown;
  okText?: string;
  /** The spaces to offer; read from the platform by default. */
  load?: () => Promise<SpaceTree>;
}

/** Answers the space saved to, or null when the dialog was cancelled. */
export function saveToSpace(options: SaveToSpaceOptions): Promise<DG.Project | null> {
  return new Promise((resolve) => {
    let saved: DG.Project | null = null;
    const dialog = new Dialog(options.title ?? 'Save to space');
    const taken = signal<{target: DG.Project, name: string} | null>(null);
    const missing = signal<string | null>(null);
    let space: SpaceInput | null = null;
    const clash = computed(() => {
      const at = taken.value;
      return at === null ? null : `${SpaceInput.describe(at.target, space?.tree.value ?? null)} already has ` +
        `"${at.name}". Rename it or pick another space.`;
    });
    let name!: TextInput;
    let description: TextArea | null = null;
    dialog.add(() => {
      const input = space = new SpaceInput({label: 'Space', value: options.space ?? null, root: options.root,
        nullable: false, load: options.load});
      name = new TextInput({label: 'Name', value: options.name ?? '', nullable: false});
      description = options.description === undefined ? null :
        new TextArea({label: 'Description', value: options.description});
      const exists = options.exists ?? spaceHasName;
      let turn = 0;
      input.effect(() => {
        const target = input.value;
        const text = name.value.trim();
        const at = ++turn;
        taken.value = null;
        if (target === null || text === '') {
          missing.value = target === null ? 'Pick a space' : 'Enter a name';
          return;
        }
        missing.value = null;
        exists(target, text).then((held) => {
          if (at === turn && held)
            taken.value = {target, name: text};
        }).catch(() => {});
      });
      return div([input, name, ...(description === null ? [] : [description]), ...(options.inputs ?? []),
        div([span(clash)], 'u2-space-clash')], 'u2-space-save');
    });
    dialog
      .okDisabled(computed(() => missing.value ?? clash.value))
      .onOK(async () => {
        const target = space!.value!;
        if (await (options.exists ?? spaceHasName)(target, name.value.trim()))
          throw new Error(`${SpaceInput.describe(target, space!.tree.peek())} already has "${name.value.trim()}".`);
        await options.save(target, {name: name.value.trim(), description: description?.value ?? ''});
        rememberSpace(target);
        saved = target;
      }, {text: options.okText ?? 'SAVE', primary: true})
      .onClose(() => resolve(saved))
      .show({width: 520});
  });
}

export interface SpaceBrowserOptions {
  /** Entity types to list (`Project`, `FuncCall`, …); everything a space holds but its subspaces by default. */
  types?: string[];
  /** Keeps only the items it accepts: the runs of one function, say. */
  filter?: (item: DG.Entity) => boolean;
  multi?: boolean;
  /** Browses this space and the spaces under it alone. */
  root?: DG.Project;
  load?: () => Promise<SpaceTree>;
  /** What a space holds itself; its children through `grok.dapi.spaces` by default. */
  items?: (space: DG.Project, types: string[]) => Promise<DG.Entity[]>;
  /** Double-click or Enter on an item. */
  onActivate?: (item: DG.Entity) => void;
}

interface Located {
  item: DG.Entity;
  space: SpaceEntry;
}

const ALL = 'all';

async function spaceItems(space: DG.Project, types: string[]): Promise<DG.Entity[]> {
  const found = await grok.dapi.spaces.id(space.id).children.ofTypes(types.join(','), true).list();
  return found.filter((e) => !(e instanceof DG.Project && e.isSpace));
}

/** The spaces and what they hold: pick a space, or "All spaces", and search what it and its subspaces hold. */
export class SpaceBrowser extends Control {
  /** The picked items, the lead first. */
  readonly selected: ReadonlySignal<DG.Entity[]>;
  readonly spaces = signal<SpaceTree | null>(null);
  private readonly _options: SpaceBrowserOptions;
  private readonly _held = new Map<string, Promise<DG.Entity[]>>();

  constructor(options: SpaceBrowserOptions = {}) {
    super();
    this._options = options;
    this.root.classList.add('u2-space-browser');
    this.root.dataset.u2 = 'space-browser';

    const scopeId = signal(options.root?.id ?? ALL);
    const query = signal('');
    const found = signal<Located[] | null>(null);
    const failed = signal<string | null>(null);

    const tree = this.build(() => new VirtualTree<SpaceEntry | null>({itemHeight: 24,
      render: (node) => divH([icon(node.data === null ? 'layer-group' : 'brackets-curly', {cls: 'u2-space-node-icon'}),
        span(node.label)], 'u2-space-node')}));
    tree.root.classList.add('u2-space-browser-tree');
    const search = this.build(() => new TextInput({inline: true, bind: query,
      placeholder: computed(() => scopeId.value === ALL ? 'Search all spaces' :
        `Search ${this.spaces.value?.get(scopeId.value)?.project.friendlyName ?? 'this space'}`)}));
    search.root.classList.add('u2-space-search');

    const renderer = handlerRenderer<DG.Entity>();
    const list = this.build(() => new VirtualList<Located>({itemHeight: 40,
      keyOf: (x) => `${x.space.project.id}/${x.item.id}`,
      labelOf: (x) => x.item.friendlyName,
      render: (x) => divH([div([renderer.icon(x.item)], 'u2-space-item-icon'),
        div([span(x.item.friendlyName, 'u2-space-item-name'),
          span(this.spaces.peek()?.path(x.space.project).join(' / ') ?? '', 'u2-space-dim')], 'u2-space-item-text'),
        span(Format.value(x.item.updatedOn ?? x.item.createdOn), 'u2-space-dim')], 'u2-space-item'),
      onEnter: (x) => this._options.onActivate?.(x.item),
    }));
    list.root.classList.add('u2-space-items');

    const shown = computed<Located[]>(() => {
      const text = query.value.trim().toLowerCase();
      return (found.value ?? []).filter((x) => text === '' || x.item.friendlyName.toLowerCase().includes(text));
    });
    list.setItems(shown);
    this.selected = computed(() => {
      const items = shown.value;
      const picked = [...list.selectedIndices.value].map((i) => items[i]?.item).filter((x): x is DG.Entity => !!x);
      return options.multi === true ? picked : picked.slice(0, 1);
    });

    const empty = div([], 'u2-space-empty-state');
    this.effect(() => {
      const items = shown.value;
      empty.hidden = found.value === null || items.length > 0;
      empty.textContent = failed.value ?? (query.value.trim() !== '' ? `Nothing matches "${query.value.trim()}".` :
        scopeId.value === ALL ? 'The spaces hold nothing to open here.' : 'This space holds nothing to open here.');
    });
    const loading = div([span('Loading…')], 'u2-space-status');
    this.effect(() => loading.hidden = found.value !== null || failed.value !== null);

    const preview = div([], 'u2-space-preview');
    this.effect(() => {
      const lead = this.selected.value[0];
      preview.replaceChildren(...(lead === undefined ? [div([span('Pick an item to see it here')], 'u2-space-dim')] :
        [entityCard(lead, {currentOnClick: false}).root]));
    });

    this.effect(() => {
      const spaces = this.spaces.value;
      if (spaces === null)
        return;
      const tops = options.root === undefined ? spaces.roots : [spaces.get(options.root.id)]
        .filter((e): e is SpaceEntry => e !== null);
      const toNode = (e: SpaceEntry): TreeNode<SpaceEntry | null> => ({id: e.project.id, label: e.project.friendlyName,
        data: e, children: e.children.length > 0 ? e.children.map(toNode) : undefined});
      const nodes = tops.map(toNode);
      tree.setRoots(options.root === undefined ? [{id: ALL, label: 'All spaces', data: null}, ...nodes] : nodes);
      void tree.expandPath([scopeId.peek()]);
    });
    tree.effect(() => {
      const node = tree.selectedNode.value;
      if (node !== null)
        scopeId.value = node.id;
    });

    this.effect(() => {
      const spaces = this.spaces.value;
      const id = scopeId.value;
      if (spaces === null)
        return;
      found.value = null;
      failed.value = null;
      list.selectedIndex.value = -1;
      const entries = id === ALL ? spaces.all() : spaces.under(id);
      Promise.all(entries.map(async (space) =>
        (await this._itemsOf(space.project)).map((item) => ({item, space}))))
        .then((lists) => {
          if (scopeId.peek() === id && !this.scope.isDisposed)
            found.value = lists.flat();
        })
        .catch((e) => {
          if (scopeId.peek() === id && !this.scope.isDisposed)
            failed.value = `Could not read the spaces: ${e?.message ?? e}`;
        });
    });

    (options.load ?? SpaceTree.load)()
      .then((spaces) => {
        if (!this.scope.isDisposed)
          this.spaces.value = spaces;
      })
      .catch((e) => failed.value = `Could not read the spaces: ${e?.message ?? e}`);

    this.root.append(div([tree.root], 'u2-space-browser-left'),
      div([search.root, loading, empty, list.root], 'u2-space-browser-middle'), preview);
  }

  // read once per space
  private _itemsOf(space: DG.Project): Promise<DG.Entity[]> {
    let held = this._held.get(space.id);
    if (held === undefined) {
      const types = this._options.types ?? [];
      const filter = this._options.filter;
      held = (this._options.items ?? spaceItems)(space, types)
        .then((items) => filter === undefined ? items : items.filter(filter));
      this._held.set(space.id, held);
    }
    return held;
  }
}

export interface OpenFromSpaceOptions extends Omit<SpaceBrowserOptions, 'onActivate'> {
  title?: string;
  /** What OPEN does with the picked items; the platform's double-click action on each by default. */
  open?: (items: DG.Entity[]) => unknown;
}

/** Answers the items opened, or an empty list when the dialog was cancelled. */
export function openFromSpace(options: OpenFromSpaceOptions = {}): Promise<DG.Entity[]> {
  return new Promise((resolve) => {
    let opened: DG.Entity[] = [];
    let browser!: SpaceBrowser;
    const open = options.open ?? ((items: DG.Entity[]) => {
      for (const item of items)
        DG.ObjectHandler.runDefaultAction(item);
    });
    const dialog = new Dialog(options.title ?? 'Open from space');
    const pick = signal<string | null>('Pick an item first');
    dialog.add(() => {
      browser = new SpaceBrowser({...options, onActivate: () => dialog.getButton('OPEN')?.click()});
      browser.effect(() => pick.value = browser.selected.value.length === 0 ? 'Pick an item first' : null);
      return browser;
    });
    dialog
      .okDisabled(pick)
      .onOK(async () => {
        opened = browser.selected.peek();
        await open(opened);
      }, {text: 'OPEN', primary: true})
      .onClose((reason) => resolve(reason === 'ok' ? opened : []))
      .show({width: 960, height: 560, resizable: true});
  });
}

