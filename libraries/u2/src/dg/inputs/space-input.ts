import type * as DG from 'datagrok-api/dg';
import {Input, InputOptions} from '../../core/input-base.js';
import {signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {div, divH, span, on} from '../../core/elements.js';
import {icon} from '../../components/display/icon.js';
import {iconButton} from '../../components/actions/buttons.js';
import {TextInput} from '../../components/inputs/text-input.js';
import {VirtualTree} from '../../components/collections/tree.js';
import type {TreeNode} from '../../components/collections/tree.js';
import {SpaceTree, SpaceEntry, createSpace, recentSpaces} from '../spaces/spaces.js';

export interface SpaceInputOptions extends InputOptions<DG.Project | null> {
  placeholder?: string;
  /** Offers this space and the spaces under it alone. */
  root?: DG.Project;
  /** The "+" buttons that create a space: at the top level and inside a space; on by default. */
  allowNew?: boolean;
  /** The spaces to offer; read from the platform by default. */
  load?: () => Promise<SpaceTree>;
}

interface Row {
  /** Null on the row of a space being named. */
  entry: SpaceEntry | null;
}

const NEW_NAME = 'New space';

export class SpaceInput extends Input<DG.Project | null, SpaceInputOptions> {
  // all assigned in createEditor, which the base constructor runs before field initializers
  private declare _tree: Signal<SpaceTree | null>;
  private declare _box: HTMLElement;
  private declare _read: Promise<SpaceTree> | null;
  private declare _close: (() => void) | null;

  constructor(options: SpaceInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'space-input';
    // the path of a given space needs the tree; a failed read is shown when the popup opens
    if (options.value)
      this.load().catch(() => {});
  }

  /** The spaces, once read. */
  get tree(): ReadonlySignal<SpaceTree | null> {
    return this._tree;
  }

  /** Reads the spaces once; a later call answers the same read. */
  load(): Promise<SpaceTree> {
    this._read ??= (this.options.load ?? SpaceTree.load)().then((tree) => {
      if (!this.scope.isDisposed)
        this._tree.value = tree;
      return tree;
    });
    return this._read;
  }

  /** The value as a path: "Chemistry / Kinase". */
  static describe(target: DG.Project | null, tree: SpaceTree | null): string {
    if (target === null)
      return '';
    return tree === null ? target.friendlyName : tree.path(target).join(' / ');
  }

  protected createEditor(): HTMLElement {
    this._tree = signal(null);
    this._read = null;
    this._close = null;
    const text = span('', 'u2-space-text');
    const clear = icon('times', {cls: 'u2-space-clear'});
    clear.title = 'Clear';
    const box = div([icon('brackets-curly', {cls: 'u2-space-icon'}), text, clear,
      icon('chevron-down', {cls: 'u2-space-caret'})], 'u2-space-box');
    box.tabIndex = 0;
    box.setAttribute('role', 'combobox');
    box.setAttribute('aria-haspopup', 'tree');
    this._box = box;
    this.scope.listen(box, 'click', (e) => {
      if ((e.target as HTMLElement).closest('.u2-space-clear'))
        this.value = null;
      else
        this._close === null ? this.open() : this._close();
    });
    this.scope.listen(box, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key === 'Enter' || key === 'ArrowDown' || key === ' ') {
        e.preventDefault();
        this.open();
      }
    });
    this.effect(() => {
      const shown = SpaceInput.describe(this.value, this.tree.value);
      text.textContent = shown === '' ? this.options.placeholder ?? 'Not in a space' : shown;
      box.classList.toggle('u2-space-empty', shown === '');
      clear.hidden = shown === '' || !this.nullable;
    });
    return box;
  }

  /** Opens the space tree under the box. */
  open(): void {
    if (this._close !== null || this.disabled)
      return;
    const scope = this.scope.child('space-popup');
    const popup = scope.build(() => this._popup());
    on(popup, OVERLAY_CLOSE_EVENT, () => {
      this._close = null;
      this._box.classList.remove('u2-space-open');
      this._box.setAttribute('aria-expanded', 'false');
      scope.dispose();
    });
    this._box.classList.add('u2-space-open');
    this._box.setAttribute('aria-expanded', 'true');
    this._close = Overlay.show(this._box, popup, {owner: this.scope, maxHeight: 420});
    queueMicrotask(() => popup.querySelector<HTMLInputElement>('.u2-space-search input')?.focus());
  }

  private _pick(target: DG.Project | null): void {
    this.value = target;
    this._close?.();
  }

  private _popup(): HTMLElement {
    const query = signal('');
    const search = new TextInput({inline: true, placeholder: 'Search spaces', bind: query});
    search.root.classList.add('u2-space-search');
    let pending: {id: string, parent: DG.Project | null} | null = null;
    let creating: Promise<void> | null = null;
    const failed = div([], 'u2-space-failed');

    const tree = new VirtualTree<Row>({
      itemHeight: 24,
      render: (node) => {
        const entry = node.data?.entry ?? null;
        const parts = [icon(entry === null ? 'plus' : 'brackets-curly', {cls: 'u2-space-node-icon'}), span(node.label)];
        if (entry !== null && entry.readOnly === null && this.options.allowNew !== false) {
          const add = iconButton('plus', () => newSpace(entry.project).catch(() => {}),
            {tooltip: `New space inside ${entry.project.friendlyName}`});
          add.classList.add('u2-space-add');
          parts.push(add);
        }
        return divH(parts, 'u2-space-node');
      },
      onRename: (node, label) => {
        if (pending !== null && node.id === pending.id)
          endNaming(label.trim());
      },
    });

    const endNaming = (name: string | null) => {
      const spaces = this.tree.peek();
      if (pending === null || spaces === null)
        return;
      const parent = pending.parent;
      pending = null;
      if (name === null || name === '') {
        rebuild(spaces);
        return;
      }
      // created right away, as a file dialog's new folder is: it stays even when the save goes elsewhere
      const before = this.value;
      creating = createSpace(name, parent).then(async (project) => {
        const entry = spaces.add(project, parent);
        // a space the user picked meanwhile stays picked, and a closed popup has no tree to show it in
        if (this.value === before)
          this.value = project;
        if (tree.scope.isDisposed)
          return;
        rebuild(spaces);
        const ids: string[] = [];
        for (let e: SpaceEntry | null = entry; e !== null; e = e.parent)
          ids.unshift(e.project.id);
        await tree.expandPath(ids);
      }).catch((e) => {
        failed.textContent = `Could not create "${name}": ${e?.message ?? e}`;
        if (!tree.scope.isDisposed)
          rebuild(spaces);
      }).finally(() => creating = null);
    };
    // the tree reports a changed name only: Enter on the name as offered keeps it, as "New folder" does
    tree.scope.listen(tree.root, 'keydown', (e) => {
      const box = (e.target as HTMLElement).closest<HTMLInputElement>('.u2-tree-rename');
      if (box !== null && (e as KeyboardEvent).key === 'Enter' && [NEW_NAME, ''].includes(box.value.trim()))
        endNaming(NEW_NAME);
    }, {capture: true});
    // a click away ends the naming with no rename; capture, as the tree's own blur takes the box out of the row
    tree.scope.listen(tree.root, 'blur', (e) => {
      const id = (e.target as HTMLElement).closest('.u2-tree-rename') ? tree.nodeForRow(e.target)?.id : undefined;
      if (id !== undefined)
        queueMicrotask(() => pending?.id === id && endNaming(null));
    }, {capture: true});
    tree.root.classList.add('u2-space-tree');

    const recents = div([], 'u2-space-recents');
    const recentHead = div([span('Recent')], 'u2-space-head');
    const allowed = (e: SpaceEntry): boolean => {
      const root = this.options.root;
      for (let c: SpaceEntry | null = e; c !== null; c = c.parent) {
        if (root === undefined || c.project.id === root.id)
          return true;
      }
      return false;
    };

    const toNode = (e: SpaceEntry, text: string, spaces: SpaceTree): TreeNode<Row> | null => {
      const kids = e.children.map((c) => toNode(c, text, spaces)).filter((n): n is TreeNode<Row> => n !== null);
      if (pending !== null && pending.parent?.id === e.project.id)
        kids.push(pendingNode());
      const hit = text === '' || spaces.path(e.project).join(' / ').toLowerCase().includes(text);
      if (!hit && kids.length === 0)
        return null;
      return {id: e.project.id, label: e.project.friendlyName, data: {entry: e}, disabled: e.readOnly,
        children: kids.length > 0 ? kids : undefined};
    };
    const pendingNode = (): TreeNode<Row> =>
      ({id: pending!.id, label: NEW_NAME, data: {entry: null}});

    const rebuild = (spaces: SpaceTree) => {
      const text = query.value.trim().toLowerCase();
      const tops = this.options.root === undefined ? spaces.roots : [spaces.get(this.options.root.id)]
        .filter((e): e is SpaceEntry => e !== null);
      const nodes = tops.map((e) => toNode(e, text, spaces)).filter((n): n is TreeNode<Row> => n !== null);
      if (pending !== null && pending.parent === null)
        nodes.push(pendingNode());
      tree.setRoots(nodes);
      const open = new Set(tree.expanded.peek());
      if (text !== '') {
        for (const e of spaces.all())
          open.add(e.project.id);
      }
      const at = this.value;
      for (let c = at === null ? null : spaces.get(at.id); c !== null; c = c.parent)
        open.add(c.project.id);
      if (pending?.parent)
        open.add(pending.parent.id);
      tree.expanded.value = open;

      const recent = recentSpaces().map((id) => spaces.get(id))
        .filter((e): e is SpaceEntry => e !== null && e.readOnly === null && allowed(e));
      recents.replaceChildren(...recent.map((e) => on(divH([icon('history', {cls: 'u2-space-node-icon'}),
        span(e.project.friendlyName), span(spaces.path(e.project).slice(0, -1).join(' / '), 'u2-space-dim')],
      'u2-space-recent'), 'click', () => this._pick(e.project))));
      recents.hidden = recentHead.hidden = text !== '' || recent.length === 0;
    };
    tree.effect(() => {
      const spaces = this.tree.value;
      void query.value;
      if (spaces !== null)
        rebuild(spaces);
    });

    const pickNode = (node: TreeNode<Row> | null) => {
      if (pending !== null || node === null || node.disabled)
        return;
      if (node.data?.entry)
        this._pick(node.data.entry.project);
    };
    on(tree.root, 'click', (e) => {
      const target = e.target as HTMLElement;
      if (!target.closest('.u2-tree-twistie, .u2-tree-rename, .u2-space-add'))
        pickNode(tree.nodeForRow(target));
    });
    on(tree.root, 'keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Enter' && !(e.target as HTMLElement).closest('.u2-tree-rename'))
        pickNode(tree.selectedNode.peek());
    });

    const status = div([span('Loading spaces…')], 'u2-space-status');
    tree.effect(() => status.hidden = this.tree.value !== null);
    this.load().catch((e) => status.textContent = `Could not read the spaces: ${e?.message ?? e}`);

    // a "+" during a create waits for it: the create's own redraw would take the new name box away
    const newSpace = async (parent: DG.Project | null) => {
      const spaces = await this.load();
      await creating;
      if (!tree.root.isConnected || pending !== null)
        return;
      pending = {id: `new:${Date.now()}`, parent};
      failed.textContent = '';
      query.value = '';
      rebuild(spaces);
      await tree.rename(pending.id);
    };

    const root = this.options.root ?? null;
    const head = divH([span(root === null ? 'Spaces' : root.friendlyName)], 'u2-space-head');
    if (this.options.allowNew !== false) {
      head.append(iconButton('plus', () => newSpace(root).catch(() => {}),
        {tooltip: root === null ? 'New top-level space' : `New space inside ${root.friendlyName}`}));
    }
    return div([search.root, recentHead, recents, head, status, failed, tree.root], 'u2-space-popup');
  }
}

export function spaceInput(label: string, options: SpaceInputOptions = {}): SpaceInput {
  return new SpaceInput({...options, label});
}
