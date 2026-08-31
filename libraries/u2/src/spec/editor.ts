/* Editing the document (DD10): every gesture — a palette drop, a grid keystroke, an AI edit — is a
   patch against the live spec. The engine validates it, mutates the node tree in place, re-renders
   the subtrees the change touched, and keeps the inverse so undo replays through the same path.
   Platform-free: the designer view and the property panel are gesture layers over this. */
import {computed, signal} from '../core/signals.js';
import type {ReadonlySignal} from '../core/signals.js';
import {SpecInstance, checkProp, htmlProps, isHtmlTag} from './spec.js';
import type {SpecEventEntry, SpecNode} from './spec.js';
import {findBindingCycle, parsePath, referencesName, referencesOf} from './path.js';
import type {ComponentMeta, SpecPropMeta} from './registry.js';

export type SpecPatch =
  | {op: 'add', parent: SpecNode, index: number, node: SpecNode}
  | {op: 'remove', node: SpecNode}
  | {op: 'move', node: SpecNode, parent: SpecNode, index: number}
  | {op: 'set-prop', node: SpecNode, name: string, value?: unknown}
  | {op: 'set-bind', node: SpecNode, name: string, path?: string}
  | {op: 'set-event', node: SpecNode, event: string, command?: SpecEventEntry}
  | {op: 'rename', node: SpecNode, name?: string}
  | {op: 'add-component', index: number, node: SpecNode}
  | {op: 'remove-component', node: SpecNode};

export interface AppliedPatch {
  /** One entry per {@link SpecEditor.apply}, several per {@link SpecEditor.applyAll}. */
  patches: SpecPatch[];
  origin: 'apply' | 'undo' | 'redo';
  /** add/remove/move/rename in ANY member — what the tree and the selection must react to. */
  structural: boolean;
}

interface UndoEntry {
  patches: SpecPatch[];
  inverses: SpecPatch[];
  origin: 'apply' | 'redo';
  /** Monotonic edit stamp: bumped on push AND on a coalescing merge, kept across undo/redo —
   * so the stack top identifies a document state, and {@link SpecEditor.dirty} is a comparison. */
  serial: number;
}

const UNDO_LIMIT = 200;
const IDENTIFIER = /^[A-Za-z_$][A-Za-z0-9_$]*$/;

export class SpecEditor {
  readonly canUndo: ReadonlySignal<boolean>;
  readonly canRedo: ReadonlySignal<boolean>;
  /** True while the document differs from the last {@link markClean} point — NOT `canUndo`:
   * undoing back to the mark is clean again, and a saved document keeps its history. */
  readonly dirty: ReadonlySignal<boolean>;
  /** The last change that went through, from any origin — the view and the panel converge on it. */
  readonly onDidApply: ReadonlySignal<AppliedPatch | null>;

  private readonly _undo: UndoEntry[] = [];
  private readonly _redo: UndoEntry[] = [];
  private readonly _version = signal(0);
  private readonly _lastApplied = signal<AppliedPatch | null>(null);
  private _serial = 0;
  private _mark = 0;

  constructor(readonly instance: SpecInstance) {
    // the stacks are plain arrays: reading the version is what makes their length reactive
    this.canUndo = computed(() => this._version.value >= 0 && this._undo.length > 0);
    this.canRedo = computed(() => this._version.value >= 0 && this._redo.length > 0);
    this.dirty = computed(() => this._version.value >= 0 && this._topSerial() !== this._mark);
    this.onDidApply = this._lastApplied;
  }

  /** Declares the current document state saved: {@link dirty} answers false until it changes. */
  markClean(): void {
    this._mark = this._topSerial();
    this._version.value++;
  }

  private _topSerial(): number {
    return this._undo.length > 0 ? this._undo[this._undo.length - 1].serial : 0;
  }

  /** Null when the patch may apply; the human-readable refusal otherwise. Unknown tags,
   * unresolvable bind paths and unknown commands pass — containing those is the renderer's job. */
  canApply(patch: SpecPatch): string | null {
    switch (patch.op) {
      case 'add': {
        const meta = this.instance.registry.get(patch.node.tag);
        if (meta?.visual === false)
          return `${patch.node.tag} is a data source — it belongs on the components tray`;
        return this._checkParent(patch.parent) ?? this._checkIndex(patch.parent, patch.index, false);
      }
      case 'add-component': {
        const components = this._components();
        if (components.includes(patch.node))
          return `"${label(patch.node)}" is already on the tray`;
        const meta = this.instance.registry.get(patch.node.tag);
        // an unknown tag is allowed on either side: containment renders it broken, and the user
        // can fix the tag in place
        if (meta && meta.visual !== false)
          return `${patch.node.tag} is a visual component — it belongs in the form`;
        return Number.isInteger(patch.index) && patch.index >= 0 && patch.index <= components.length ?
          null : `index ${patch.index} is out of range for the tray (0..${components.length})`;
      }
      case 'remove-component':
        return this._components().includes(patch.node) ? null :
          `"${label(patch.node)}" is not on the tray`;
      case 'remove':
        return this._checkStructural(patch.node);
      case 'move': {
        const refusal = this._checkStructural(patch.node) ?? this._checkParent(patch.parent) ??
          this._checkIndex(patch.parent, patch.index, this.instance.parentOf(patch.node) === patch.parent);
        if (refusal)
          return refusal;
        return patch.node === patch.parent || SpecInstance.contains(patch.node, patch.parent) ?
          `"${label(patch.node)}" cannot move into itself` : null;
      }
      case 'set-prop': {
        const found = this._prop(patch.node, patch.name);
        if (!found)
          return `${patch.node.tag} has no prop "${patch.name}"`;
        if (!this._rendered(found.own ? patch.node : this.instance.parentOf(patch.node)!))
          return `"${label(patch.node)}" is not rendered`;
        return patch.value === undefined ? null : checkProp(found.prop, patch.value);
      }
      case 'set-bind': {
        if (!this._rendered(patch.node))
          return `"${label(patch.node)}" is not rendered`;
        if (patch.path === undefined)
          return null;
        const dot = patch.name.indexOf('.');
        if (dot >= 0) {
          const head = patch.name.slice(0, dot);
          const found = this._prop(patch.node, head);
          if (!found)
            return `${patch.node.tag} has no prop "${head}"`;
          if (!found.prop.subBindable)
            return `prop "${head}" does not take sub-binds`;
        } else {
          const found = this._prop(patch.node, patch.name);
          if (!found)
            return `${patch.node.tag} has no prop "${patch.name}"`;
          if (found.prop.subBindable && !found.prop.bindable)
            return `prop "${patch.name}" takes sub-binds only`;
        }
        if (!patch.path.startsWith('$.'))
          return `a bind path must start with "$." — got "${patch.path}"`;
        if (parsePath(patch.path) === null)
          return `"${patch.path}" is not a valid bind path`;
        return this._checkCycle(patch.node, patch.name, patch.path);
      }
      case 'set-event':
        if (!this._rendered(patch.node))
          return `"${label(patch.node)}" is not rendered`;
        return patch.command === undefined ? null : SpecEditor._checkEvent(patch.command);
      case 'rename': {
        if (!this._rendered(patch.node))
          return `"${label(patch.node)}" is not rendered`;
        if (patch.name === undefined)
          return null;
        if (!IDENTIFIER.test(patch.name))
          return `"${patch.name}" is not a valid name: letters, digits, _ and $, not starting with a digit`;
        const owner = this._named().get(patch.name);
        if (owner !== undefined && owner !== patch.node)
          return `"${patch.name}" is already taken`;
        // the rewrite below moves every `$.<name>`/`cmd:<name>.` reference onto the new name, and
        // those namespaces are shared with the context: renaming onto a data key or a command would
        // capture references that were never this node's, and undo could not tell them apart
        return SpecEditor._inContext(this.instance, patch.name) ?
          `"${patch.name}" is already a context data key or command` : null;
      }
    }
  }

  /** Throws what {@link canApply} refuses. */
  apply(patch: SpecPatch): void {
    const refusal = this.canApply(patch);
    if (refusal)
      throw new Error(`u2 editor: ${refusal}`);
    const inverse = this._inverse(patch);
    // history first: a re-render that throws mid-way still leaves the stack able to unwind
    const top = this._undo[this._undo.length - 1];
    if (top && top.origin === 'apply' && top.patches.length === 1 &&
      SpecEditor._mergeable(top.patches[0], patch)) {
      top.patches = [patch];
      // the merged entry now lands on a different document state than the one it may be marked at
      top.serial = ++this._serial;
    } else {
      this._undo.push({patches: [patch], inverses: [inverse], origin: 'apply', serial: ++this._serial});
      if (this._undo.length > UNDO_LIMIT)
        this._undo.shift();
    }
    this._redo.length = 0;
    this._execute(patch);
    this._publish([patch], 'apply');
  }

  /** Validates and applies in order as ONE undo entry; a mid-sequence refusal unwinds the
   * already-applied members (inverses, reverse order) and throws. Never coalesces. */
  applyAll(patches: SpecPatch[]): void {
    if (patches.length === 0)
      return;
    const applied: SpecPatch[] = [];
    const inverses: SpecPatch[] = [];
    for (const patch of patches) {
      const refusal = this.canApply(patch);
      if (refusal) {
        this._unwind(inverses);
        throw new Error(`u2 editor: ${refusal}`);
      }
      inverses.push(this._inverse(patch));
      try {
        this._execute(patch);
      } catch (e) {
        // the throwing member's mutation landed before its re-render — its inverse unwinds it too
        this._unwind(inverses);
        throw e;
      }
      applied.push(patch);
    }
    this._undo.push({patches: applied, inverses, origin: 'apply', serial: ++this._serial});
    if (this._undo.length > UNDO_LIMIT)
      this._undo.shift();
    this._redo.length = 0;
    this._publish(applied, 'apply');
  }

  /** Reverse order: the inverse of each member was computed against the members before it. */
  private _unwind(inverses: SpecPatch[]): void {
    for (let i = inverses.length - 1; i >= 0; i--)
      this._execute(inverses[i]);
  }

  undo(): boolean {
    const entry = this._undo.pop();
    if (entry === undefined)
      return false;
    // reverse order: the inverse of each member was computed against the members before it
    const executed: SpecPatch[] = [];
    for (let i = entry.inverses.length - 1; i >= 0; i--) {
      this._execute(entry.inverses[i]);
      executed.push(entry.inverses[i]);
    }
    this._redo.push(entry);
    this._publish(executed, 'undo');
    return true;
  }

  redo(): boolean {
    const entry = this._redo.pop();
    if (entry === undefined)
      return false;
    for (const patch of entry.patches)
      this._execute(patch);
    this._undo.push({...entry, origin: 'redo'});
    this._publish(entry.patches, 'redo');
    return true;
  }

  /** Serializable address forms for tools outside the object graph (AI patches, P6). */
  nodePath(node: SpecNode): number[] | null {
    const path: number[] = [];
    const walk = (at: SpecNode): boolean => {
      if (at === node)
        return true;
      const children = at.children ?? [];
      for (let i = 0; i < children.length; i++) {
        path.push(i);
        if (walk(children[i]))
          return true;
        path.pop();
      }
      return false;
    };
    return walk(this.instance.spec.root) ? path : null;
  }

  nodeAt(path: number[]): SpecNode | null {
    let at: SpecNode = this.instance.spec.root;
    for (const index of path) {
      const child = at.children?.[index];
      if (child === undefined)
        return null;
      at = child;
    }
    return at;
  }

  /** `textInput1`-style: `base` plus the smallest counter free across the whole spec. */
  uniqueName(base: string): string {
    return this.uniqueNames()(base);
  }

  /** A run of names that also reserves what it hands out: the nodes of a seeded subtree are not in
   * the document yet, so {@link uniqueName} alone would hand every sibling the same one. */
  uniqueNames(): (base: string) => string {
    const taken = new Set(this._named().keys());
    return (base) => {
      const name = SpecEditor._free(base, taken);
      taken.add(name);
      return name;
    };
  }

  /** Deep copy with every name re-uniquified; not inserted — pair it with an `add` patch. */
  duplicate(node: SpecNode): SpecNode {
    const clone = JSON.parse(JSON.stringify(node)) as SpecNode;
    const taken = this._named();
    const rename = (at: SpecNode): void => {
      if (at.name !== undefined) {
        at.name = SpecEditor._free(at.name.replace(/\d+$/, ''), taken);
        taken.set(at.name, at);
      }
      for (const child of at.children ?? [])
        rename(child);
    };
    rename(clone);
    return clone;
  }

  /** Whatever references the patched node captured its signals at its own render, and a rebuilt —
   * or removed — node leaves them wired to the corpse: the references render again after it, in
   * the same batch. Read after the mutation, so a rename finds them under the new name. */
  private _execute(patch: SpecPatch): void {
    const targets = this._targets(patch);
    targets.push(...this._mutate(patch));
    if (patch.node.name !== undefined)
      targets.push(...referencesOf(this.instance.spec, patch.node.name));
    this.instance.rerenderAll(targets);
  }

  private _components(): SpecNode[] {
    return this.instance.spec.components ?? [];
  }

  /** What the change is anchored to, read before the mutation: the parent for anything structural
   * or for a prop the parent declares, the node itself otherwise. */
  private _targets(patch: SpecPatch): SpecNode[] {
    switch (patch.op) {
      case 'add':
        return [patch.parent];
      // the tray has no visual anchor: the mutation mounts and releases the component itself
      case 'add-component':
      case 'remove-component':
        return [];
      case 'remove':
        return [this.instance.parentOf(patch.node)!];
      case 'move':
        return [this.instance.parentOf(patch.node)!, patch.parent];
      case 'set-prop':
        return [this._prop(patch.node, patch.name)!.own ? patch.node : this.instance.parentOf(patch.node)!];
      default:
        return [patch.node];
    }
  }

  private _publish(patches: SpecPatch[], origin: 'apply' | 'undo' | 'redo'): void {
    this._version.value++;
    this._lastApplied.value = {patches, origin, structural: patches.some((patch) =>
      patch.op !== 'set-prop' && patch.op !== 'set-bind' && patch.op !== 'set-event')};
  }

  private _inverse(patch: SpecPatch): SpecPatch {
    switch (patch.op) {
      case 'add':
        return {op: 'remove', node: patch.node};
      case 'add-component':
        return {op: 'remove-component', node: patch.node};
      case 'remove-component':
        return {op: 'add-component', index: this._components().indexOf(patch.node), node: patch.node};
      case 'remove':
      case 'move': {
        const parent = this.instance.parentOf(patch.node)!;
        const index = (parent.children ?? []).indexOf(patch.node);
        return patch.op === 'remove' ?
          {op: 'add', parent, index, node: patch.node} :
          {op: 'move', node: patch.node, parent, index};
      }
      case 'set-prop':
        return {op: 'set-prop', node: patch.node, name: patch.name, value: patch.node.props?.[patch.name]};
      case 'set-bind':
        return {op: 'set-bind', node: patch.node, name: patch.name, path: patch.node.bind?.[patch.name]};
      case 'set-event':
        return {op: 'set-event', node: patch.node, event: patch.event, command: patch.node.on?.[patch.event]};
      case 'rename':
        return {op: 'rename', node: patch.node, name: patch.node.name};
    }
  }

  /** Unlike an unresolvable path, a cycle can livelock re-runs — a validation concern, not
   * containment. The would-be spec is probed in place: the entry is set, the whole document
   * checked, the entry restored exactly. Only a loop through this node is this patch's doing —
   * a pre-existing loop in a hand-written spec is the renderer's containment, not a refusal. */
  private _checkCycle(node: SpecNode, name: string, path: string): string | null {
    if (node.name === undefined)
      return null;
    const had = Object.prototype.hasOwnProperty.call(node.bind ?? {}, name);
    const old = node.bind?.[name];
    SpecEditor._entry(node, 'bind', name, path);
    const loop = findBindingCycle(this.instance.spec);
    SpecEditor._entry(node, 'bind', name, had ? old : undefined);
    return loop !== null && loop.includes(node.name) ?
      `binding "${path}" would create a cycle: ${loop.join(' → ')}` : null;
  }

  private static _checkEvent(entry: SpecEventEntry): string | null {
    const cmd = typeof entry === 'string' ? entry : entry?.cmd;
    if (typeof cmd !== 'string' || !cmd.startsWith('cmd:')) {
      // no angle brackets: the platform balloon renders HTML, and its sanitizer would eat
      // a "<name>" placeholder together with everything after it
      return 'an event must name a command as \'cmd:\' followed by the command name — got ' +
        `'${typeof cmd === 'string' ? cmd : JSON.stringify(entry)}'`;
    }
    if (typeof entry !== 'string' && entry.args !== undefined &&
      (typeof entry.args !== 'object' || entry.args === null || Array.isArray(entry.args)))
      return 'event args must be a plain object of argument values';
    return null;
  }

  private _checkStructural(node: SpecNode): string | null {
    if (node === this.instance.spec.root)
      return 'the root node cannot be removed or moved';
    if (this._components().includes(node))
      return `"${label(node)}" is a component — the tray has its own add and remove`;
    const parent = this.instance.parentOf(node);
    if (parent === null)
      return `"${label(node)}" is not in the spec`;
    return this._rendered(parent) ? null : `"${label(parent)}" is not rendered`;
  }

  private _checkParent(parent: SpecNode): string | null {
    const meta = this.instance.registry.get(parent.tag);
    const accepts = meta ? meta.acceptsChildren === true || meta.adopt !== undefined ||
      meta.createWithChildren !== undefined : isHtmlTag(parent.tag);
    if (!accepts)
      return `${parent.tag} takes no children`;
    return this._rendered(parent) ? null : `"${label(parent)}" is not rendered`;
  }

  private _checkIndex(parent: SpecNode, index: number, moving: boolean): string | null {
    const last = (parent.children ?? []).length - (moving ? 1 : 0);
    if (Number.isInteger(index) && index >= 0 && index <= last)
      return null;
    return `index ${index} is out of range for ${parent.tag} (0..${last})`;
  }

  /** The prop the renderer would validate this name against: the component's own, then — as
   * `_prop` does — one the parent declares for its children. */
  private _prop(node: SpecNode, name: string): {prop: SpecPropMeta, own: boolean} | null {
    const meta = this.instance.registry.get(node.tag);
    const own = meta ? meta.props.find((p) => p.name === name) :
      isHtmlTag(node.tag) ? htmlProps(node.tag).find((p) => p.name === name) : undefined;
    if (own)
      return {prop: own, own: true};
    const parent = this.instance.parentOf(node);
    const childProp = parent ?
      this.instance.registry.get(parent.tag)?.childProps?.find((p) => p.name === name) : undefined;
    return childProp ? {prop: childProp, own: false} : null;
  }

  private _rendered(node: SpecNode): boolean {
    return this.instance.nodes().has(node);
  }

  private static _inContext(instance: SpecInstance, name: string): boolean {
    const owns = Object.prototype.hasOwnProperty;
    return owns.call(instance.ctx.data, name) || owns.call(instance.ctx.commands, name);
  }

  private _named(): Map<string, SpecNode> {
    const names = new Map<string, SpecNode>();
    const walk = (node: SpecNode): void => {
      if (node.name !== undefined && !names.has(node.name))
        names.set(node.name, node);
      for (const child of node.children ?? [])
        walk(child);
    };
    // one namespace across both trees: a bind path names a component and a node the same way
    walk(this.instance.spec.root);
    for (const component of this._components())
      walk(component);
    return names;
  }

  private static _free(base: string, taken: {has(name: string): boolean}): string {
    for (let i = 1; ; i++) {
      const name = `${base}${i}`;
      if (!taken.has(name))
        return name;
    }
  }

  /** Answers the render targets the mutation itself discovers: a rename rewrites the references
   * other nodes hold, and those nodes are wired to what they referenced — they must come back too. */
  private _mutate(patch: SpecPatch): SpecNode[] {
    switch (patch.op) {
      case 'add':
        SpecEditor._insert(patch.parent, patch.index, patch.node);
        return [];
      case 'add-component': {
        const spec = this.instance.spec;
        if (spec.components)
          spec.components.splice(patch.index, 0, patch.node);
        else
          spec.components = [patch.node];
        this.instance.mountComponent(patch.node);
        return [];
      }
      case 'remove-component': {
        const components = this.instance.spec.components!;
        this.instance.unmountComponent(patch.node);
        components.splice(components.indexOf(patch.node), 1);
        // like an emptied `children`, an empty tray goes with its last entry so `dump()` stays exact
        if (components.length === 0)
          delete this.instance.spec.components;
        return [];
      }
      case 'remove':
        SpecEditor._detach(this.instance.parentOf(patch.node)!, patch.node);
        return [];
      case 'move':
        SpecEditor._detach(this.instance.parentOf(patch.node)!, patch.node);
        SpecEditor._insert(patch.parent, patch.index, patch.node);
        return [];
      case 'set-prop':
        SpecEditor._entry(patch.node, 'props', patch.name, patch.value);
        return [];
      case 'set-bind':
        SpecEditor._entry(patch.node, 'bind', patch.name, patch.path);
        return [];
      case 'set-event':
        SpecEditor._entry(patch.node, 'on', patch.event, patch.command);
        return [];
      case 'rename': {
        const old = patch.node.name;
        const touched: SpecNode[] = [];
        if (old !== undefined && patch.name !== undefined) {
          touched.push(...renameRefs(this.instance.spec.root, old, patch.name));
          for (const component of this.instance.spec.components ?? [])
            touched.push(...renameRefs(component, old, patch.name));
        }
        if (patch.name === undefined)
          delete patch.node.name;
        else
          patch.node.name = patch.name;
        return touched.filter((node) => this._rendered(node));
      }
    }
  }

  private static _insert(parent: SpecNode, index: number, node: SpecNode): void {
    if (parent.children)
      parent.children.splice(index, 0, node);
    else
      parent.children = [node];
  }

  private static _detach(parent: SpecNode, node: SpecNode): void {
    const children = parent.children!;
    children.splice(children.indexOf(node), 1);
    if (children.length === 0)
      delete parent.children;
  }

  /** An emptied container goes with its last entry, so `dump()` never grows `props: {}` noise. */
  private static _entry(node: SpecNode, key: 'props' | 'bind' | 'on', name: string, value: unknown): void {
    const container = node[key] as Record<string, unknown> | undefined;
    if (value !== undefined) {
      if (container)
        container[name] = value;
      else
        (node as unknown as Record<string, unknown>)[key] = {[name]: value};
      return;
    }
    if (!container)
      return;
    delete container[name];
    if (Object.keys(container).length === 0)
      delete node[key];
  }

  /** Typing into one field is one undo step: the same op on the same key merges, keeping the
   * original inverse and the latest forward patch. */
  private static _mergeable(a: SpecPatch, b: SpecPatch): boolean {
    if (a.op !== b.op || a.node !== b.node)
      return false;
    if (a.op === 'set-prop' && b.op === 'set-prop')
      return a.name === b.name;
    if (a.op === 'set-bind' && b.op === 'set-bind')
      return a.name === b.name;
    if (a.op === 'set-event' && b.op === 'set-event')
      return a.event === b.event;
    return false;
  }
}

/** camelCased tag suffix: 'u2-text-input' → 'textInput'; an HTML tag answers itself. */
export function nameForTag(tag: string): string {
  const parts = tag.startsWith('u2-') ? tag.slice(3).split('-') : [tag];
  return parts[0] + parts.slice(1).map((part) => part.charAt(0).toUpperCase() + part.slice(1)).join('');
}

/** The node a designer inserts for `tag`: the meta's `defaults` deep-cloned, plus — when the meta
 * declares a `label` prop the defaults leave unset — a label humanized from the unique name, so a
 * dropped input has a visible identity from the first frame (the Delphi convention). A multi-host
 * (tabs, accordion) seeds its `defaultChildren` the same way, each one named through `unique`: they
 * ride the container's own `add` patch, so one undo takes the whole insert back. `unique` names
 * them — a tag that declares no `defaultChildren` (every non-visual source among them) needs none. */
export function seedNode(meta: ComponentMeta | undefined, tag: string, name: string,
  unique?: (base: string) => string): SpecNode {
  const node: SpecNode = {tag, name};
  const props: Record<string, unknown> =
    meta?.defaults ? JSON.parse(JSON.stringify(meta.defaults)) : {};
  if (meta?.props.some((p) => p.name === 'label') && props.label === undefined)
    props.label = humanize(name);
  if (Object.keys(props).length > 0)
    node.props = props;
  if (meta?.defaultChildren && unique) {
    node.children = meta.defaultChildren.map((child) =>
      nameAll(JSON.parse(JSON.stringify(child)) as SpecNode, unique));
  }
  return node;
}

/** A seeded child is addressable exactly like an inserted one — in the tree, by `data-u2-name` and
 * by the rename gesture — so every node of the cloned subtree gets a free name off its own tag. */
function nameAll(node: SpecNode, unique: (base: string) => string): SpecNode {
  node.name = unique(nameForTag(node.tag));
  for (const child of node.children ?? [])
    nameAll(child, unique);
  return node;
}

/** Sentence case off the camelCased unique name: 'textInput1' → 'Text input 1'. */
function humanize(name: string): string {
  const words = name.replace(/([A-Z])/g, ' $1').replace(/(\d+)$/, ' $1').toLowerCase();
  return words.charAt(0).toUpperCase() + words.slice(1);
}

/** First-segment rewrite of `$.old…` bind paths and `cmd:old.…` event values — the references a
 * component name owns (DD5). A bare `cmd:old` names a context command, not a node, and is left
 * alone. Answers the nodes whose references moved: what they resolve against changed under them,
 * so they have to be rendered again. */
export function renameRefs(root: SpecNode, oldName: string, newName: string): SpecNode[] {
  const touched: SpecNode[] = [];
  const renamedCmd = (command: string): string | null =>
    command.startsWith(`cmd:${oldName}.`) ?
      `cmd:${newName}.${command.slice(oldName.length + 5)}` : null;
  const walk = (node: SpecNode): void => {
    let moved = false;
    for (const [name, path] of Object.entries(node.bind ?? {})) {
      const next = renamedPath(path, oldName, newName);
      if (next !== null) {
        node.bind![name] = next;
        moved = true;
      }
    }
    for (const [event, entry] of Object.entries(node.on ?? {})) {
      if (typeof entry === 'string') {
        const next = renamedCmd(entry);
        if (next !== null) {
          node.on![event] = next;
          moved = true;
        }
        continue;
      }
      const cmd = renamedCmd(entry.cmd);
      if (cmd !== null) {
        entry.cmd = cmd;
        moved = true;
      }
      for (const [key, value] of Object.entries(entry.args ?? {})) {
        const next = typeof value === 'string' ? renamedPath(value, oldName, newName) : null;
        if (next !== null) {
          entry.args![key] = next;
          moved = true;
        }
      }
    }
    if (moved)
      touched.push(node);
    for (const child of node.children ?? [])
      walk(child);
  };
  walk(root);
  return touched;
}

function renamedPath(path: string, oldName: string, newName: string): string | null {
  return referencesName(path, oldName) ? `$.${newName}${path.slice(oldName.length + 2)}` : null;
}

function label(node: SpecNode): string {
  return node.name ?? node.tag;
}
