/* The editing verbs over the current selection (D-8): one list feeding the ribbon buttons and every
   context menu, each enabled exactly where the engine would take it. */
import * as grok from 'datagrok-api/grok';
import {Component} from '../../core/component.js';
import {Reason} from '../../core/reason.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {actionsMenu} from '../../components/actions/actions.js';
import type {Action} from '../../components/actions/actions.js';
import type {TreeNode} from '../../components/collections/tree.js';
import {SpecInstance} from '../../spec/spec.js';
import type {SpecNode} from '../../spec/spec.js';
import type {DesignerAction} from '../../spec/registry.js';
import {nameForTag, seedNode} from '../../spec/editor.js';
import type {SpecEditor, SpecPatch} from '../../spec/editor.js';
import {refreshSource} from './source-status.js';

/** The editing verbs, in the order the ribbon shows them: `actions()` names them from here, and
 * the ribbon reads their enablement back by position. */
export const ACTIONS = ['Delete', 'Move Up', 'Move Down', 'Duplicate'] as const;
export const [DELETE, MOVE_UP, MOVE_DOWN, DUPLICATE] = ACTIONS;
export const ACTION_ICONS: Record<string, string> = {
  [DELETE]: 'trash-alt', [MOVE_UP]: 'arrow-up', [MOVE_DOWN]: 'arrow-down', [DUPLICATE]: 'clone'};

export interface ActionsHost {
  readonly editor: ReadonlySignal<SpecEditor | undefined>;
  instance(): SpecInstance | undefined;
  selected(): SpecNode | null;
  multi(): SpecNode[];
  select(node: SpecNode): void;
  pending(node: SpecNode | null): void;
  apply(patch: SpecPatch): boolean;
  refocus(): void;
}

export class DesignerActions {
  /** Why a reorder or a duplicate waits for a single selection. */
  static readonly ONE_AT_A_TIME = 'Select one component for this';
  /** Why the form's root or a data source on the tray is not reordered or duplicated. */
  static readonly NO_PARENT = 'It has no parent to move or copy it within';

  /** The verbs that may ask the user first: they refocus when they are done, not when they start. */
  private readonly _deferred = new WeakSet<Action>();

  constructor(private readonly _host: ActionsHost) {}

  /** A multi-selection deletes as one compound patch; the reorder verbs wait for a single node (F5). */
  actions(): Action[] {
    const node = this._host.selected();
    const editor = this._host.editor.peek();
    const instance = this._host.instance();
    if (!node || !editor || !instance)
      return [];
    if (this._host.multi().length > 1) {
      const cover = this.cover();
      return [
        {text: DELETE, icon: ACTION_ICONS[DELETE],
          disabled: cover.map((member) => editor.canApply({op: 'remove', node: member}))
            .find((r) => r !== null) ?? null,
          run: () => this._deleteMulti(cover)},
        {text: MOVE_UP, icon: ACTION_ICONS[MOVE_UP], disabled: DesignerActions.ONE_AT_A_TIME, run: () => {}},
        {text: MOVE_DOWN, icon: ACTION_ICONS[MOVE_DOWN], disabled: DesignerActions.ONE_AT_A_TIME, run: () => {}},
        {text: DUPLICATE, icon: ACTION_ICONS[DUPLICATE], disabled: DesignerActions.ONE_AT_A_TIME, run: () => {}},
      ];
    }
    const parent = instance.parentOf(node);
    const siblings = parent?.children ?? [];
    const index = siblings.indexOf(node);
    // a tray component has no parent, so the reorder verbs disable themselves; only Delete needs
    // to know which of the two removals it is
    return [
      {text: DELETE, icon: ACTION_ICONS[DELETE], disabled: editor.canApply(this._removalOf(node)),
        run: () => this._delete(node)},
      {text: MOVE_UP, icon: ACTION_ICONS[MOVE_UP],
        disabled: parent === null ? DesignerActions.NO_PARENT : index > 0 ? null : 'Already the first',
        run: () => void this._host.apply({op: 'move', node, parent: parent!, index: index - 1})},
      {text: MOVE_DOWN, icon: ACTION_ICONS[MOVE_DOWN],
        disabled: parent === null ? DesignerActions.NO_PARENT : index < siblings.length - 1 ? null : 'Already the last',
        run: () => void this._host.apply({op: 'move', node, parent: parent!, index: index + 1})},
      {text: DUPLICATE, icon: ACTION_ICONS[DUPLICATE], disabled: parent === null ? DesignerActions.NO_PARENT : null,
        run: () => this._duplicate(node)},
      // after the fixed four: the ribbon reads only those by position, so these are menu-only
      ...this._refreshAction(node),
      ...(instance.registry.get(node.tag)?.designerActions ?? []).map((action): Action => {
        const verb: Action = {text: action.text, icon: action.icon,
          run: () => void this._designerAction(node, action)};
        this._deferred.add(verb);
        return verb;
      }),
    ];
  }

  /** Running a source once is an explicit user act — the one thing design time never does by
   * itself (DD9). Offered on every component that declares the function. */
  private _refreshAction(node: SpecNode): Action[] {
    return this._source(node) === undefined ? [] :
      [{text: 'Refresh', icon: 'sync', run: () => this.refresh(node)}];
  }

  /** The run design time makes only when asked (DD9): the chip's Refresh, and a dropped file —
   * the drop IS the ask, and a source with no visible pulse reads as broken. */
  refresh(node: SpecNode): void {
    const source = this._source(node);
    if (source)
      void refreshSource(source, node.name ?? node.tag);
  }

  private _source(node: SpecNode): Component | undefined {
    const built = this._host.instance()?.nodes().get(node);
    return this._isComponent(node) && Component.is(built) &&
      built.getFunctions().some((f) => f.name === 'refresh') ? built : undefined;
  }

  private _isComponent(node: SpecNode): boolean {
    return this._host.instance()?.spec.components?.includes(node) === true;
  }

  private _removalOf(node: SpecNode): SpecPatch {
    return this._isComponent(node) ? {op: 'remove-component', node} : {op: 'remove', node};
  }

  /** A meta action describes a change against the node's current state (registry.ts); the view
   * turns it into a patch — an added child is seeded and named exactly like a palette drop. */
  private async _designerAction(node: SpecNode, action: DesignerAction): Promise<void> {
    const instance = this._host.instance();
    if (!instance)
      return;
    try {
      const produced = await action.produce(node, instance.nodes().get(node));
      const editor = this._host.editor.peek();
      if (produced === null || !editor)
        return;
      if (produced.op === 'set-prop') {
        this._host.apply({op: 'set-prop', node, name: produced.name, value: produced.value});
        return;
      }
      const unique = editor.uniqueNames();
      const tag = produced.node.tag;
      const seeded = seedNode(instance.registry.get(tag), tag, unique(nameForTag(tag)), unique);
      // the action's own props win over the seeded ones — the pane title it counted out is the point
      const child: SpecNode = {...seeded, ...produced.node, name: seeded.name};
      const props = {...seeded.props, ...produced.node.props};
      if (Object.keys(props).length > 0)
        child.props = props;
      this._host.pending(child);
      this._host.apply({op: 'add', parent: node, index: (node.children ?? []).length, node: child});
    } catch (e) {
      grok.shell.warning(`${action.text}: ${e instanceof Error ? e.message : String(e)}`);
    } finally {
      this._host.refocus();
    }
  }

  run(name: string): void {
    const action = this.actions().find((a) => a.text === name);
    if (action && Reason.now(action.disabled) === null) {
      action.run();
      if (!this._deferred.has(action))
        this._host.refocus();
    }
  }

  handBack(actions: Action[]): Action[] {
    return actions.map((a) => this._deferred.has(a) ? a : {...a, run: () => {
      a.run();
      this._host.refocus();
    }});
  }

  private _delete(node: SpecNode): void {
    const instance = this._host.instance()!;
    // a removed component leaves nothing beside it to select: the form's root takes over
    this._host.pending(this._isComponent(node) ? instance.spec.root : instance.parentOf(node));
    this._host.apply(this._removalOf(node));
  }

  /** The minimal cover of the multi-selection: a member whose ancestor is also selected goes with
   * that ancestor's remove, so it must not get one of its own. */
  cover(): SpecNode[] {
    const multi = this._host.multi();
    return multi.filter((node) =>
      !multi.some((other) => other !== node && SpecInstance.contains(other, node)));
  }

  /** One compound patch, one undo entry — each remove pre-checked through the same warn funnel a
   * single delete uses (sound: the cover's removes are disjoint). */
  private _deleteMulti(cover: SpecNode[]): void {
    const editor = this._host.editor.peek();
    const instance = this._host.instance();
    if (!editor || !instance)
      return;
    for (const node of cover) {
      const refusal = editor.canApply({op: 'remove', node});
      if (refusal) {
        grok.shell.warning(refusal);
        return;
      }
    }
    // the lead's own parent may be going too — the survivor is outside what the cover removes
    const lead = this._host.selected()!;
    const anchor = cover.find((node) => node === lead || SpecInstance.contains(node, lead))!;
    this._host.pending(instance.parentOf(anchor));
    editor.applyAll(cover.map((node): SpecPatch => ({op: 'remove', node})));
  }

  private _duplicate(node: SpecNode): void {
    const parent = this._host.instance()!.parentOf(node)!;
    const clone = this._host.editor.peek()!.duplicate(node);
    const index = (parent.children ?? []).indexOf(node) + 1;
    this._host.pending(clone);
    this._host.apply({op: 'add', parent, index, node: clone});
  }

  rowActions(row: TreeNode<SpecNode>): Action[] {
    const multi = this._host.multi();
    // a right-click inside the multi-selection keeps it — the menu acts on the whole set
    if (row.data && !(multi.length > 1 && multi.includes(row.data)))
      this._host.select(row.data);
    return this.handBack(this.actions());
  }

  /** The design surface never shows the browser menu, nor the platform's — even when it has no
   * menu of its own to offer. */
  onContextMenu(e: MouseEvent, node: SpecNode | null): void {
    e.preventDefault();
    e.stopPropagation();
    const multi = this._host.multi();
    if (node && !(multi.length > 1 && multi.includes(node)))
      this._host.select(node);
    const actions = this.actions();
    if (actions.length === 0)
      return;
    actionsMenu(this.handBack(actions)).show({x: e.clientX, y: e.clientY});
  }
}
