/* dg-ui/1: JSON in, a live component tree out, and back through dump(). A malformed node becomes
   a visible placeholder — the rest of the tree always renders. Specs never carry code: events
   name context commands, and that is the only executable form. */
import {signal, Signal, isWritableSignal} from '../core/signals.js';
import {bindTypeOf} from '../core/widget-like.js';
import {Scope} from '../core/scope.js';
import {Component, Control} from '../core/component.js';
import {ComponentMeta, ComponentStart, SpecPropMeta, Registry, registry as globalRegistry,
  SPEC_SCHEMA} from './registry.js';
import {findBindingCycle, parsePath, referencesOf} from './path.js';
import {BindProp, BindSource, isBindSource} from './bind-source.js';
import {APPEARANCE_PROPS, applyAppearance} from './appearance.js';

/** An event's wiring: `'cmd:name'`, or the structured form carrying arguments — values that are
 * strings starting `$.` are paths resolved when the event fires, `$$.` escapes a literal `$.`,
 * everything else is a literal. Still JSON, still code-free. */
export type SpecEventEntry = string | {cmd: string, args?: Record<string, unknown>};

export interface SpecNode {
  tag: string;
  /** Unique within the spec: what selection, patches, code-behind and automation address. */
  name?: string;
  props?: Record<string, unknown>;
  /** `{propName: '$.path'}` — binds a prop to what the path resolves to. */
  bind?: Record<string, string>;
  /** `{eventName: 'cmd:name'}` — the only executable form a spec has. */
  on?: Record<string, SpecEventEntry>;
  children?: SpecNode[];
}

export interface Spec {
  $schema: string;
  root: SpecNode;
  /** The non-visual tray (Q1): data sources and state, same node shape, no children. Built before
   * the visual tree — a bind resolves against them — and started after it. */
  components?: SpecNode[];
}

export interface SpecInstanceOptions {
  /** The designer's mode: sources build for preview, never for effect (DD9). */
  designTime?: boolean;
}

export interface SpecContextOptions {
  data?: Record<string, unknown>;
  commands?: Record<string, () => void>;
  /** The platform-function tier of `cmd:` (Q6); absent — gallery, headless — that tier resolves
   * nothing and the handler warns. */
  callFunction?: (name: string, args: Record<string, unknown>) => Promise<unknown>;
}

/** What one node rendered: the scope everything under it hangs on, and the handle that disposes it
 * and drops the disposer its owner keeps — so re-rendering the same node many times leaves no dead
 * closures on the parent. */
interface Mount {
  scope: Scope;
  release(): void;
}

/** A bind whose target is a visual node declared later in the document. `link`: the component took
 * a proxy signal at construction, and the link to the real signal waits for the render pass to end.
 * `wire`: an HTML node's own live prop — {@link SpecInstance._flushLinks} resolves the target and
 * runs {@link Deferred.wire}, so the element is never re-rendered. `rerender`: a re-render-tier
 * bind built with `undefined` — the flush re-renders the node with the resolved value baked in
 * ({@link SpecInstance._rebaked}), then wires the follow once the rebuilt value converges. */
interface Deferred {
  node: SpecNode;
  /** Identity guard, as before: a re-rendered node's stale entries are skipped. */
  built: Component | HTMLElement;
  kind: 'link' | 'wire' | 'rerender';
  name: string;
  path: string;
  twoWay: boolean;
  /** `wire` only: installs the in-place follower over the resolved source. */
  wire?: (source: Signal<unknown>) => void;
}

// module-level is safe only because builds are synchronous and every builder deletes its entry first
const childBindValues = new WeakMap<SpecNode, Record<string, unknown>>();

/** The value a parent reads off a child's declared prop (a pane's `title`): the resolved bind
 * when the child carries one, the literal otherwise. */
export function childProp(node: SpecNode | undefined, name: string): unknown {
  const bound = node ? childBindValues.get(node)?.[name] : undefined;
  return bound !== undefined ? bound : node?.props?.[name];
}

const HTML_TAGS = new Set(['div', 'span', 'p', 'h1', 'h2', 'h3', 'a', 'img']);
const HTML_PROPS = new Set(['text', 'cls', 'href', 'src']);
const JSON_TYPES = new Set(['object', 'string', 'number', 'boolean']);

/** What `$.path` and `cmd:name` resolve against. Plain data values are wrapped into signals,
 * so a spec binds to them and sees live updates either way. */
export class SpecContext {
  readonly data: Record<string, Signal<unknown> | BindSource> = {};
  readonly commands: Record<string, () => void>;
  readonly callFunction?: (name: string, args: Record<string, unknown>) => Promise<unknown>;

  constructor(options: SpecContextOptions = {}) {
    const data = options.data ?? {};
    for (const key of Object.keys(data)) {
      const value = data[key];
      this.data[key] = value instanceof Signal ? value as Signal<unknown> :
        isBindSource(value) ? value : signal(value);
    }
    this.commands = options.commands ?? {};
    this.callFunction = options.callFunction;
  }

  /** What the binding picker roots at the context: one entry per data key, the type inferred from
   * the current value, sources walkable. Allocates nothing. */
  bindProps(): BindProp[] {
    return Object.entries(this.data).map(([name, value]) => value instanceof Signal ?
      {name, type: bindTypeOf(value.peek()), writable: isWritableSignal(value)} :
      {name, walkable: true});
  }
}

/** A rendered spec: a scope per node — so one subtree can be released and rendered again alone —
 * plus the node → built map that makes {@link dump} read live state and {@link nodeAt} hit-test
 * the DOM. */
export class SpecInstance extends Control {
  readonly ctx: SpecContext;

  private readonly _spec: Spec;
  private readonly _registry: Registry;
  private readonly _nodes = new Map<SpecNode, Component | HTMLElement>();
  private readonly _named = new Map<string, Component | HTMLElement>();
  private readonly _elements = new Map<Element, SpecNode>();
  private readonly _mounts = new Map<SpecNode, Mount>();
  private readonly _warned = new Set<string>();
  /** Components a re-render has mounted but not started yet: the same two phases construction has,
   * for the same reason — a param bound to an input must resolve against the input the re-render
   * leaves behind, not the one it replaced. */
  private readonly _pending: SpecNode[] = [];
  private readonly _deferred: Deferred[] = [];
  /** Per node, the resolved value a forward-ref re-render bind is rebuilt with. Under an ADOPTING
   * render root the flush's re-render re-creates the target too, so the rebuilt node still sees a
   * forward ref — building with the baked value and wiring at the next flush is what terminates
   * that loop with the value resolved. `rounds` caps re-bakes for a target whose value is a fresh
   * object on every build, which `Object.is` can never call converged. */
  private readonly _rebaked = new WeakMap<SpecNode, Map<string, {value: unknown, rounds: number}>>();
  private readonly _rerenderQueue = new Set<SpecNode>();
  private _rerenderScheduled = false;
  private _cycle: string[] | null;
  private _designTime: boolean;
  private _batching = false;

  constructor(spec: Spec, ctx: SpecContext, reg: Registry = globalRegistry,
    options: SpecInstanceOptions = {}) {
    super();
    this._spec = spec;
    this.ctx = ctx;
    this._registry = reg;
    this._designTime = options.designTime ?? false;
    this._cycle = findBindingCycle(spec);
    this.root.classList.add('u2-spec');
    // two phases (Q1): the tray is built first, so a visual bind resolves against it, and started
    // last, so a source param can bind to an input declared anywhere in the form
    this.runInScope(() => {
      for (const node of spec.components ?? [])
        this._mountComponent(node);
      this.root.append(SpecInstance._element(this._render(spec.root)));
      this._flushLinks();
      for (const node of spec.components ?? [])
        this._start(node);
    });
  }

  dispose(): void {
    super.dispose();
    this._nodes.clear();
    this._named.clear();
    this._elements.clear();
    this._mounts.clear();
    this._deferred.length = 0;
  }

  get spec(): Spec {
    return this._spec;
  }

  get registry(): Registry {
    return this._registry;
  }

  /** The designer's Design/Run toggle (DD9). Not a patch: the mode is view state, so the edit
   * history is untouched — but every tray component is built again for the new mode, and
   * {@link rerenderAll}'s dependent expansion brings everything bound to one back with it
   * (a visual node captured the old signals), before the sources start. */
  setDesignTime(x: boolean): void {
    if (x === this._designTime)
      return;
    this._designTime = x;
    this.rerenderAll(this._spec.components ?? []);
  }

  get designTime(): boolean {
    return this._designTime;
  }

  /** The document as it is: bound props keep their `bind`, plain ones their literal — a Run-mode
   * edit of a control is never folded in; the designer's Design-mode panel edits are patches and
   * therefore in the document. Rendering the result reproduces this UI. */
  dump(): object {
    const components = this._spec.components ?? [];
    const out: Record<string, unknown> = {$schema: SPEC_SCHEMA};
    // authoring order: what the form binds to is declared before the form
    if (components.length > 0)
      out.components = components.map((node) => this._dump(node));
    out.root = this._dump(this._spec.root);
    return out;
  }

  /** Every node the spec built, tray components, plain HTML and placeholders alike. */
  nodes(): ReadonlyMap<SpecNode, Component | HTMLElement> {
    return this._nodes;
  }

  node(name: string): Component | HTMLElement | undefined {
    return this._named.get(name);
  }

  /** The element a build occupies — none for a tray component, which never reaches the DOM. */
  static elementOf(built: Component | HTMLElement): HTMLElement | undefined {
    return Control.is(built) ? built.root : Component.is(built) ? undefined : built;
  }

  /** The nearest spec node owning `el` — the designer's hit-test. */
  nodeAt(el: Element): SpecNode | null {
    for (let at: Element | null = el; at; at = at.parentElement) {
      const node = this._elements.get(at);
      if (node)
        return node;
    }
    return null;
  }

  /** Walks from the spec root; null for the root or an unreachable node. */
  parentOf(node: SpecNode): SpecNode | null {
    return node === this._spec.root ? null : SpecInstance._parentOf(this._spec.root, node);
  }

  /** Whether `descendant` sits anywhere under `node` — the cycle check a move and a drop need. */
  static contains(node: SpecNode, descendant: SpecNode): boolean {
    return SpecInstance._parentOf(node, descendant) !== null;
  }

  /** Full binding resolution (Q2). First segment: a named component → a named visual node → a
   * context data entry; then a {@link BindSource} walk, `bindStep('')` — the default binding —
   * at path exhaustion. Throws with a precise message; render-time containment turns that into
   * the node's placeholder. */
  resolveBinding(path: string): {signal: Signal<unknown>, writable: boolean} {
    const segments = parsePath(path);
    if (segments === null)
      throw new Error(`"${path}" is not a valid bind path`);
    let at = this._firstStep(segments[0], path);
    for (const segment of segments.slice(1)) {
      if (at instanceof Signal)
        throw new Error(`"${path}": cannot walk into "${segment}" — the path continues past a signal`);
      const next = at.bindStep(segment);
      if (next === null)
        throw new Error(`"${path}": no binding step "${segment}"`);
      at = next;
    }
    // depth-bounded, like the picker's own walk: a default step that answers its own source would
    // otherwise spin here and hang the renderer
    for (let depth = 3; !(at instanceof Signal); depth--) {
      const source: BindSource = at;
      const fallback = depth > 0 ? source.bindStep('') : null;
      if (fallback === null) {
        const choices = source.bindProps().map((p) => p.name).join(', ');
        throw new Error(`"${path}" needs one more step — one of: ${choices}`);
      }
      at = fallback;
    }
    return {signal: at, writable: isWritableSignal(at)};
  }

  private _firstStep(first: string, path: string): Signal<unknown> | BindSource {
    const component = this._spec.components?.find((c) => c.name === first);
    if (component) {
      const built = this._nodes.get(component);
      if (built !== undefined && isBindSource(built))
        return built;
      throw new Error(`component "${first}" is not built`);
    }
    const named = this._named.get(first);
    if (named !== undefined) {
      if (isBindSource(named))
        return named;
      if (SpecInstance._isPlaceholder(named))
        throw new Error(`component "${first}" is declared but not built`);
      throw new Error(`"${first}" is not a bind source`);
    }
    const data = this.ctx.data[first];
    if (data !== undefined)
      return data;
    if (SpecInstance._declared(this._spec.root, first))
      throw new Error(`component "${first}" is declared but not built`);
    throw new Error(`nothing bound at "${path}"`);
  }

  private static _eachName(node: SpecNode, fn: (name: string) => void): void {
    if (node.name !== undefined)
      fn(node.name);
    for (const child of node.children ?? [])
      SpecInstance._eachName(child, fn);
  }

  /** Whether the document declares a node of this name at all, built or not. */
  private static _declared(root: SpecNode, name: string): boolean {
    return root.name === name || (root.children ?? []).some((child) => SpecInstance._declared(child, name));
  }

  /** A later-declared visual node: the one first segment {@link resolveBinding} cannot answer yet. */
  private _forwardRef(path: string): boolean {
    const first = parsePath(path)?.[0];
    return first !== undefined && !this._named.has(first) && this.ctx.data[first] === undefined &&
      !(this._spec.components ?? []).some((c) => c.name === first) && SpecInstance._declared(this._spec.root, first);
  }

  /** The render pass is over: every forward reference links to what the document built — or
   * warns, leaving the referencing node built and unlinked (the target's own placeholder names
   * the real problem). */
  private _flushLinks(): void {
    const rerenders: SpecNode[] = [];
    for (const d of this._deferred.splice(0)) {
      const mount = this._mounts.get(d.node);
      if (mount === undefined || mount.scope.isDisposed)
        continue;
      if (this._nodes.get(d.node) !== d.built)
        continue;
      try {
        if (d.kind === 'wire') {
          const {signal: source} = this.resolveBinding(d.path);
          Scope.runWith(mount.scope, () => d.wire!(source));
        } else if (d.kind === 'rerender') {
          const {signal: source} = this.resolveBinding(d.path);
          const baked = this._rebaked.get(d.node);
          const entry = baked?.get(d.name);
          // converged on the target's current value — or round-capped: wire the follow, done
          if (entry !== undefined && (Object.is(entry.value, source.value) || entry.rounds >= 2)) {
            baked!.delete(d.name);
            Scope.runWith(mount.scope, () => this._watchRerender(d.node, [[d.name, source]]));
          } else {
            const map = baked ?? new Map();
            this._rebaked.set(d.node, map);
            map.set(d.name, {value: source.value, rounds: (entry?.rounds ?? 0) + 1});
            rerenders.push(d.node);
          }
        } else {
          const {signal: source, writable} = this.resolveBinding(d.path);
          if (d.twoWay && !writable)
            this._warn(`bind "${d.path}" is read-only — edits will not flow back`);
          Scope.runWith(mount.scope, () => (d.built as Control).link(d.name, source, d.twoWay && writable));
        }
      } catch (e) {
        this._rebaked.get(d.node)?.delete(d.name);
        this._warn(`${d.node.tag}: ${e instanceof Error ? e.message : String(e)}`);
      }
    }
    if (rerenders.length > 0)
      this.rerenderAll(rerenders);
  }

  /** Disposes what `node` rendered — its whole subtree, scopes, components, listeners and map
   * entries — and renders it anew in place. The node must currently be rendered. */
  rerender(node: SpecNode): void {
    const old = this._nodes.get(node);
    if (old === undefined)
      throw new Error('u2 spec: rerender of a node that is not rendered');
    this._cycle = findBindingCycle(this._spec);
    if (this._isComponent(node)) {
      // nothing to replace: a tray component holds no place in the DOM
      this._unmount(node);
      this._mountComponent(node);
      // once per flush, whatever rebuilt it: `start()` is not idempotent — a second one registers
      // the source's effects again and doubles its server traffic for the session
      if (!this._pending.includes(node))
        this._pending.push(node);
      if (!this._batching)
        this._flushStarts();
      return;
    }
    const oldEl = SpecInstance.elementOf(old)!;
    const parent = this.parentOf(node);
    const carried = this._carried(node);
    this._unmount(node);
    const scope = parent ? this._mounts.get(parent)!.scope : this.scope;
    const meta = parent ? this._registry.get(parent.tag) : undefined;
    const built = Scope.runWith(scope, () => this._render(node, meta));
    for (const [at, carry, prev] of carried) {
      const next = this._nodes.get(at);
      if (next !== undefined)
        carry(at, prev, next);
    }
    oldEl.replaceWith(SpecInstance._element(built));
    if (!this._batching)
      this._flushStarts();
  }

  /** View state the tags under `node` want moved across a rebuild ({@link ComponentMeta.carry}),
   * captured before the unmount that would lose it. */
  private _carried(node: SpecNode): [SpecNode, NonNullable<ComponentMeta['carry']>,
    Component | HTMLElement][] {
    const out: [SpecNode, NonNullable<ComponentMeta['carry']>, Component | HTMLElement][] = [];
    const walk = (at: SpecNode): void => {
      const carry = this._registry.get(at.tag)?.carry;
      const built = this._nodes.get(at);
      if (carry !== undefined && built !== undefined)
        out.push([at, carry, built]);
      for (const child of at.children ?? [])
        walk(child);
    };
    walk(node);
    return out;
  }

  /** Renders each node again at its render root, nested targets collapsed into the outermost —
   * what a change touching several nodes at once has to do. Nodes that are not rendered are
   * skipped: a reference may point at a subtree the same change has just dropped. */
  rerenderAll(nodes: SpecNode[]): void {
    const roots: SpecNode[] = [];
    const add = (node: SpecNode): void => {
      if (!this._nodes.has(node))
        return;
      const root = this.renderRootOf(node);
      if (!roots.includes(root))
        roots.push(root);
    };
    for (const node of nodes)
      add(node);
    // a rebuilt root re-creates every named node under it — promotion included — and whatever
    // references one from OUTSIDE stays wired to the corpse: the dependents render too, after
    // their source's root (fixed point; each name expands once, and names are finite)
    const expanded = new Set<string>();
    for (let i = 0; i < roots.length; i++) {
      SpecInstance._eachName(roots[i], (name) => {
        if (expanded.has(name))
          return;
        expanded.add(name);
        for (const dependent of referencesOf(this._spec, name))
          add(dependent);
      });
    }
    this._batching = true;
    try {
      const tops = roots.filter((root, i) =>
        !roots.some((other, j) => j !== i && SpecInstance.contains(other, root)));
      // tray sources before visual roots: a visual node rebuilt first would capture the signals
      // of a source the same batch is about to re-mount
      for (const root of tops.sort((a, b) => Number(this._isComponent(b)) - Number(this._isComponent(a))))
        this.rerender(root);
    } finally {
      this._batching = false;
    }
    // the closing act of every patch: what it mounted starts against the tree it ends up with
    this._flushStarts();
  }

  /** Promotion (D-3): a parent that wires its children into itself — `Form.add`, a tab strip
   * consuming them at construction — must be rebuilt whole; replacing one child element under it
   * would corrupt its internals. A tray component has no parent and answers itself. */
  renderRootOf(node: SpecNode): SpecNode {
    let at = node;
    for (;;) {
      const parent = this.parentOf(at);
      const meta = parent ? this._registry.get(parent.tag) : undefined;
      if (!parent || !meta || (meta.adopt === undefined && meta.createWithChildren === undefined))
        return at;
      at = parent;
    }
  }

  /** Builds a tray component into a rendered instance — the engine's `add-component`. Started by
   * the `rerenderAll` that closes the patch, once whatever references it has rendered again. */
  mountComponent(node: SpecNode): void {
    this._cycle = findBindingCycle(this._spec);
    this._mountComponent(node);
    this._pending.push(node);
  }

  /** Releases a tray component and everything it owns — the engine's `remove-component`. */
  unmountComponent(node: SpecNode): void {
    this._unmount(node);
  }

  private _isComponent(node: SpecNode): boolean {
    return this._spec.components?.includes(node) === true;
  }

  /** The tray's mount: the same scope bookkeeping a visual node gets, minus the DOM half — no
   * element, nothing appended (Q1). A build that throws maps the standard placeholder, so a broken
   * source is counted, selectable and undoable like any other node. */
  private _mountComponent(node: SpecNode): void {
    // owned by the instance, never by whatever scope happens to be ambient: a patch committed
    // while the context panel builds its form would hand a live data source to the panel's scope,
    // and the next `disposePanel()` would take it down
    const scope = this._mount(node, this.scope);
    Scope.runWith(scope, () => {
      const built = this._buildComponent(node);
      this._nodes.set(node, built);
      if (node.name === undefined)
        this._warn(`${node.tag}: a component without a name cannot be bound to`);
      else
        this._name(node, built);
    });
  }

  private _flushStarts(): void {
    // links before starts: a source param bound to an input that is itself forward-bound sees a live value
    this._flushLinks();
    const pending = this._pending.splice(0);
    for (const node of pending) {
      // a component the same patch removed again has nothing to start
      if (this._nodes.has(node))
        this._start(node);
    }
  }

  /** Phase two (Q1): the tray is built and the form is up — a source kicks its auto-run and
   * resolves the params bound to inputs anywhere in it. */
  private _start(node: SpecNode): void {
    const built = this._nodes.get(node) as Partial<ComponentStart> | undefined;
    if (typeof built?.start !== 'function')
      return;
    try {
      built.start();
    } catch (e) {
      this._warn(`${node.tag}: ${e instanceof Error ? e.message : String(e)}`);
    }
  }

  private _buildComponent(node: SpecNode): Component | HTMLElement {
    try {
      const meta = this._registry.get(node.tag);
      if (!meta?.createComponent) {
        throw new Error(meta ? 'is a visual component — it belongs in the form' :
          'is not a registered component');
      }
      this._children(node, meta, false);
      this._noCycle(node);
      const props = this._props(node, meta);
      const deferredRerenders: [string, string][] = [];
      const rerenderBinds: [string, Signal<unknown>][] = [];
      for (const [name, path] of Object.entries(node.bind ?? {})) {
        const dot = name.indexOf('.');
        if (dot >= 0) {
          if (!SpecInstance._prop(meta, name.slice(0, dot))?.subBindable)
            throw new Error(`prop "${name.slice(0, dot)}" does not take sub-binds`);
          continue; // delivered raw through env.subBinds, as today
        }
        const prop = SpecInstance._prop(meta, name);
        if (!prop)
          throw new Error(`has no prop "${name}" to bind`);
        if (prop.bindable || prop.subBindable)
          continue; // the source resolves these itself at start() (filter, id)
        const rerendered = this._rerenderTier(node, name, path, deferredRerenders, rerenderBinds);
        if (rerendered !== undefined)
          props[name] = rerendered.value;
      }
      const built = meta.createComponent(props, {
        designTime: this._designTime,
        subBinds: {...node.bind},
        resolve: (path) => {
          try {
            return this.resolveBinding(path).signal;
          } catch (e) {
            this._warn(`${node.tag}: ${e instanceof Error ? e.message : String(e)}`);
            return null;
          }
        },
      });
      for (const [name, path] of deferredRerenders)
        this._deferred.push({node, built, kind: 'rerender', name, path, twoWay: false});
      this._watchRerender(node, rerenderBinds);
      return built;
    } catch (e) {
      return SpecInstance._placeholder(`${node.tag}: ${e instanceof Error ? e.message : String(e)}`);
    }
  }

  private _render(node: SpecNode, parent?: ComponentMeta): Control | HTMLElement {
    const scope = this._mount(node);
    return Scope.runWith(scope, () => {
      const built = this._build(node, parent);
      const el = SpecInstance._element(built);
      this._nodes.set(node, built);
      this._elements.set(el, node);
      if (node.name !== undefined) {
        // the automation id: what a test, a tutorial or an agent addresses the node by
        el.dataset.u2Name = node.name;
        this._name(node, built);
      }
      return built;
    });
  }

  /** The scope everything under one node hangs on, and the release that also drops the disposer
   * its owner keeps — so re-rendering the same node many times leaves no dead closures behind. */
  private _mount(node: SpecNode, owner: Scope = Scope.ambient ?? this.scope): Scope {
    const scope = new Scope();
    const disposer = () => scope.dispose();
    owner.own(disposer);
    this._mounts.set(node, {scope, release: () => {
      scope.dispose();
      owner.disown(disposer);
    }});
    return scope;
  }

  private _name(node: SpecNode, built: Component | HTMLElement): void {
    if (this._named.has(node.name!))
      this._warn(`${node.tag}: duplicate name "${node.name}"`);
    else
      this._named.set(node.name!, built);
  }

  /** Releasing one node's scope cascades through every scope below it, so the maps are cleaned by
   * the scopes that died — not by walking `children`, which a patch may already have spliced. */
  private _unmount(node: SpecNode): void {
    this._mounts.get(node)!.release();
    const gone = new Set<Component | HTMLElement>();
    for (const [at, mount] of Array.from(this._mounts)) {
      if (!mount.scope.isDisposed)
        continue;
      const built = this._nodes.get(at);
      if (built !== undefined) {
        gone.add(built);
        this._nodes.delete(at);
        const el = SpecInstance.elementOf(built);
        if (el)
          this._elements.delete(el);
      }
      this._mounts.delete(at);
    }
    // by identity, not by name: a rename patch has already moved the name off the node
    for (const [name, built] of Array.from(this._named)) {
      if (gone.has(built))
        this._named.delete(name);
    }
  }

  /** `parent` is the meta of the component this node hangs under: a prop the parent declares in
   * `childProps` — a pane title — validates against it instead of counting as unknown. */
  private _build(node: SpecNode, parent?: ComponentMeta): Control | HTMLElement {
    try {
      const commands = SpecInstance._commands(node);
      const meta = this._registry.get(node.tag);
      const built = meta ? this._component(node, meta, parent) : this._html(node, parent);
      for (const [event, command] of commands)
        this._listen(node, built, event, command);
      return built;
    } catch (e) {
      return SpecInstance._placeholder(`${node.tag}: ${e instanceof Error ? e.message : String(e)}`);
    }
  }

  private _component(node: SpecNode, meta: ComponentMeta, parent?: ComponentMeta): Control {
    if (meta.visual === false)
      throw new Error('is a data source — it belongs on the components tray');
    this._noCycle(node);
    const props = this._props(node, meta, parent);

    childBindValues.delete(node);
    const twoWay: [string, Signal<unknown>, boolean][] = [];
    const deferredLinks: [string, string, boolean][] = [];
    const deferredRerenders: [string, string][] = [];
    const rerenderBinds: [string, Signal<unknown>][] = [];
    for (const [name, path] of Object.entries(node.bind ?? {})) {
      const dot = name.indexOf('.');
      // a sub-bind ('params.days') delivers raw through the components engine's env, never here
      if (dot >= 0 && SpecInstance._prop(meta, name.slice(0, dot))?.subBindable)
        continue;
      const prop = SpecInstance._prop(meta, name, parent);
      if (!prop)
        throw new Error(`has no prop "${name}" to bind`);
      if (prop.bindable) {
        if (this._forwardRef(path)) {
          props[name] = signal(undefined);
          deferredLinks.push([name, path, prop.twoWay === true]);
          continue;
        }
        const {signal: source, writable} = this.resolveBinding(path);
        props[name] = source;
        if (prop.twoWay) {
          if (!writable)
            this._warn(`bind "${path}" is read-only — edits will not flow back`);
          twoWay.push([name, source, writable]);
        }
      } else {
        const rerendered = this._rerenderTier(node, name, path, deferredRerenders, rerenderBinds);
        if (rerendered !== undefined) {
          props[name] = rerendered.value;
          if (SpecInstance._childProp(parent, name) === prop)
            childBindValues.set(node, {...childBindValues.get(node), [name]: rerendered.value});
        }
      }
    }

    const component = meta.createWithChildren ?
      meta.createWithChildren(props, this._children(node, meta, true), node.children ?? []) :
      this._adopt(node, meta, meta.create!(props));
    for (const [name, source, writable] of twoWay)
      component.link(name, source, writable);
    for (const [name, path, twoWay] of deferredLinks)
      this._deferred.push({node, built: component, kind: 'link', name, path, twoWay});
    for (const [name, path] of deferredRerenders)
      this._deferred.push({node, built: component, kind: 'rerender', name, path, twoWay: false});
    this._watchRerender(node, rerenderBinds);
    return component;
  }

  /** One scope-owned effect per re-render-bound node: it only ENQUEUES — the rerender that
   * disposes this very effect runs in a microtask, after the signal batch and the effect callback
   * have completed. The re-render registers a fresh effect over the freshly resolved sources. */
  private _watchRerender(node: SpecNode, binds: [string, Signal<unknown>][]): void {
    if (binds.length === 0)
      return;
    let first = true;
    Scope.ambient!.effect(() => {
      for (const [, source] of binds)
        void source.value;
      if (first) {
        first = false;
        return;
      }
      this._queueRerender(node);
    });
  }

  private _queueRerender(node: SpecNode): void {
    this._rerenderQueue.add(node);
    if (this._rerenderScheduled)
      return;
    this._rerenderScheduled = true;
    queueMicrotask(() => {
      this._rerenderScheduled = false;
      const queued = Array.from(this._rerenderQueue);
      this._rerenderQueue.clear();
      if (this.scope.isDisposed)
        return;
      this.rerenderAll(queued);
    });
  }

  /** The re-render-tier tail every non-live bind takes: a forward ref defers to the flush and
   * answers its bake, if a prior round left one; a resolvable path drops the bake, resolves and
   * watches. Undefined = a forward ref with no bake yet — nothing to build with. */
  private _rerenderTier(node: SpecNode, name: string, path: string,
    deferred: [string, string][], binds: [string, Signal<unknown>][]): {value: unknown} | undefined {
    if (this._forwardRef(path)) {
      deferred.push([name, path]);
      const entry = this._rebaked.get(node)?.get(name);
      return entry === undefined ? undefined : {value: entry.value};
    }
    this._rebaked.get(node)?.delete(name);
    const {signal: source} = this.resolveBinding(path);
    binds.push([name, source]);
    return {value: source.value};
  }

  /** The node's props, each validated against what the tag declares — an unknown one is a warning
   * and passes through, a wrong type throws and the node becomes a placeholder. */
  private _props(node: SpecNode, meta: ComponentMeta, parent?: ComponentMeta): Record<string, unknown> {
    const props: Record<string, unknown> = {};
    for (const [name, value] of Object.entries(node.props ?? {})) {
      const prop = SpecInstance._prop(meta, name, parent);
      if (!prop)
        this._warn(`${node.tag}: unknown prop "${name}"`);
      props[name] = prop ? SpecInstance._checked(prop, value) : value;
    }
    return props;
  }

  private _noCycle(node: SpecNode): void {
    if (node.name !== undefined && this._cycle !== null && this._cycle.includes(node.name))
      throw new Error(`binding cycle: ${this._cycle.join(' → ')}`);
  }

  /** Fully live (UB-4): every own prop a plain element declares — text, cls, href/src and the
   * appearance group — is bound through an engine-owned effect writing the DOM in place, so the
   * element the node built is the element it keeps. */
  private _html(node: SpecNode, parent?: ComponentMeta): HTMLElement {
    if (!HTML_TAGS.has(node.tag))
      throw new Error('neither a registered component nor a supported HTML tag');
    const el = document.createElement(node.tag);
    let prevCls: string[] = [];
    const applyHtml = (name: string, value: unknown): void => {
      if (name === 'text')
        el.textContent = value == null ? '' : String(value);
      else if (name === 'cls') {
        el.classList.remove(...prevCls);
        prevCls = String(value ?? '').split(/\s+/).filter((c) => c);
        el.classList.add(...prevCls);
      } else if (value == null || value === '')
        el.removeAttribute(name);
      else
        el.setAttribute(name, String(value));
    };
    const appearance: Record<string, unknown> = {};
    for (const [name, value] of Object.entries(node.props ?? {})) {
      if (!HTML_PROPS.has(name)) {
        const prop = APPEARANCE_PROPS.find((p) => p.name === name);
        if (prop) {
          appearance[name] = SpecInstance._checked(prop, value);
          continue;
        }
        if (!SpecInstance._childProp(parent, name))
          this._warn(`${node.tag}: unknown prop "${name}"`);
        continue;
      }
      if (typeof value !== 'string')
        throw new Error(`prop "${name}" expects string`);
      applyHtml(name, value);
    }
    childBindValues.delete(node);
    const scope = Scope.ambient!;
    // literals first: a bound appearance prop's effect installs below and wins over a same-named one
    applyAppearance(el, appearance, APPEARANCE_PROPS, scope);
    const ownProps = htmlProps(node.tag);
    const rerenderBinds: [string, Signal<unknown>][] = [];
    const deferredRerenders: [string, string][] = [];
    const follow = (name: string, source: Signal<unknown>): void => {
      if (APPEARANCE_PROPS.some((p) => p.name === name))
        applyAppearance(el, {[name]: source}, APPEARANCE_PROPS, scope);
      else
        scope.effect(() => applyHtml(name, source.value));
    };
    for (const [name, path] of Object.entries(node.bind ?? {})) {
      const own = ownProps.some((p) => p.name === name);
      if (!own && !SpecInstance._childProp(parent, name))
        throw new Error(`has no prop "${name}" to bind`);
      if (!own) {
        const rerendered = this._rerenderTier(node, name, path, deferredRerenders, rerenderBinds);
        if (rerendered !== undefined)
          childBindValues.set(node, {...childBindValues.get(node), [name]: rerendered.value});
        continue;
      }
      if (this._forwardRef(path)) {
        // wired in place at the flush — an element's own props are live, never re-rendered
        this._deferred.push({node, built: el, kind: 'wire', name, path, twoWay: false,
          wire: (source) => follow(name, source)});
        continue;
      }
      follow(name, this.resolveBinding(path).signal);
    }
    for (const [name, path] of deferredRerenders)
      this._deferred.push({node, built: el, kind: 'rerender', name, path, twoWay: false});
    this._watchRerender(node, rerenderBinds);
    for (const child of this._children(node, undefined, true))
      el.append(SpecInstance._element(child));
    return el;
  }

  private _adopt(node: SpecNode, meta: ComponentMeta, component: Control): Control {
    const accepts = meta.acceptsChildren === true || meta.adopt !== undefined;
    const children = this._children(node, meta, accepts);
    for (let i = 0; i < children.length; i++) {
      if (meta.adopt)
        meta.adopt(component, children[i], i);
      else
        component.root.append(SpecInstance._element(children[i]));
    }
    return component;
  }

  private _children(node: SpecNode, meta: ComponentMeta | undefined, accepts: boolean):
      (Control | HTMLElement)[] {
    const children = node.children ?? [];
    if (children.length > 0 && !accepts) {
      this._warn(`${node.tag}: takes no children — ${children.length} dropped`);
      return [];
    }
    return children.map((child) => this._render(child, meta));
  }

  /** A component node listens through its own event surface — component events and, on a u2
   * control, the DOM events of that name; a plain HTML node through the DOM alone. */
  private _listen(node: SpecNode, built: Control | HTMLElement, event: string, entry: SpecEventEntry): void {
    const handler = () => this._fire(node, entry);
    const scope = Scope.ambient ?? this.scope;
    if (Control.is(built)) {
      const subscription = built.onEvent(event).subscribe(handler);
      scope.own(() => subscription.unsubscribe());
      return;
    }
    built.addEventListener(event, handler);
    scope.own(() => built.removeEventListener(event, handler));
  }

  /** The `cmd:` tiers (Q6), classified by syntax at fire time — commands and functions may appear
   * after render. `#` is reserved for code-behind; `:` names a platform function; `.` names a
   * function of a named component; a bare name is a context command, then a platform function. */
  private _fire(node: SpecNode, entry: SpecEventEntry): void {
    const name = (typeof entry === 'string' ? entry : entry.cmd).slice(4);
    let args: Record<string, unknown>;
    try {
      args = typeof entry === 'string' ? {} : this._args(entry.args);
    } catch (e) {
      this._warn(`${node.tag}: ${e instanceof Error ? e.message : String(e)}`);
      return;
    }
    if (name.startsWith('#')) {
      this._warn(`${node.tag}: "cmd:${name}" is reserved — code-behind commands arrive with dg-form`);
      return;
    }
    if (name.includes(':')) {
      if (this.ctx.callFunction)
        this.ctx.callFunction(name, args);
      else
        this._warn(`${node.tag}: no function runner for "${name}"`);
      return;
    }
    if (name.includes('.')) {
      const dot = name.indexOf('.');
      const owner = name.slice(0, dot);
      const fn = name.slice(dot + 1);
      const component = this._named.get(owner);
      const found = Component.is(component) ?
        component.getFunctions().find((f) => f.name === fn) : undefined;
      if (found)
        found.apply(args);
      else
        this._warn(`${node.tag}: no function "${fn}" on "${owner}"`);
      return;
    }
    const run = this.ctx.commands[name];
    if (typeof run === 'function')
      run();
    else if (this.ctx.callFunction)
      this.ctx.callFunction(name, args);
    else
      this._warn(`${node.tag}: no command "${name}"`);
  }

  /** Argument values, resolved when the event fires: `$.` paths peeked — untracked, a handler
   * must not become a dependency — `$$.` unescaped once, everything else as-is. */
  private _args(args: Record<string, unknown> | undefined): Record<string, unknown> {
    const out: Record<string, unknown> = {};
    for (const [key, value] of Object.entries(args ?? {})) {
      out[key] = typeof value !== 'string' ? value :
        value.startsWith('$.') ? this.resolveBinding(value).signal.peek() :
          value.startsWith('$$.') ? value.slice(1) : value;
    }
    return out;
  }

  private _dump(node: SpecNode): object {
    const out: Record<string, unknown> = {tag: node.tag};
    if (node.name !== undefined)
      out.name = node.name;
    if (node.props)
      out.props = {...node.props};
    if (node.bind)
      out.bind = {...node.bind};
    if (node.on) {
      const on: Record<string, SpecEventEntry> = {};
      for (const [event, entry] of Object.entries(node.on))
        on[event] = typeof entry === 'string' ? entry : JSON.parse(JSON.stringify(entry));
      out.on = on;
    }
    if (node.children)
      out.children = node.children.map((child) => this._dump(child));
    return out;
  }

  private _warn(message: string): void {
    if (this._warned.has(message))
      return;
    this._warned.add(message);
    console.warn(`u2 spec: ${message}`);
  }

  private static _commands(node: SpecNode): [string, SpecEventEntry][] {
    const commands: [string, SpecEventEntry][] = [];
    for (const [event, entry] of Object.entries(node.on ?? {})) {
      const cmd = typeof entry === 'string' ? entry : entry?.cmd;
      if (typeof cmd !== 'string' || !cmd.startsWith('cmd:'))
        throw new Error(`on.${event} must be 'cmd:' followed by a command name — a spec never carries code`);
      commands.push([event, entry]);
    }
    return commands;
  }

  private static _prop(meta: ComponentMeta, name: string, parent?: ComponentMeta): SpecPropMeta | undefined {
    return meta.props.find((p) => p.name === name) ?? SpecInstance._childProp(parent, name);
  }

  private static _childProp(parent: ComponentMeta | undefined, name: string): SpecPropMeta | undefined {
    return parent?.childProps?.find((p) => p.name === name);
  }

  private static _element(built: Control | HTMLElement): HTMLElement {
    return Control.is(built) ? built.root : built;
  }

  private static _checked(prop: SpecPropMeta, value: unknown): unknown {
    const error = checkProp(prop, value);
    if (error)
      throw new Error(error);
    return value;
  }

  private static _parentOf(at: SpecNode, node: SpecNode): SpecNode | null {
    for (const child of at.children ?? []) {
      if (child === node)
        return at;
      const found = SpecInstance._parentOf(child, node);
      if (found)
        return found;
    }
    return null;
  }

  private static _placeholder(message: string): HTMLElement {
    const el = document.createElement('div');
    el.className = 'u2-spec-error';
    el.textContent = message;
    return el;
  }

  private static _isPlaceholder(built: Component | HTMLElement): boolean {
    return !Component.is(built) && built.classList.contains('u2-spec-error');
  }
}

const warnedTags = new Set<string>();

/** One line per distinct message for the life of the page — the dedupe {@link SpecInstance} does
 * per instance, for what is built outside one: a tray component is constructed again on every
 * design/run toggle and every prop edit, and its complaint is about the document, not the mount. */
export function specWarn(message: string): void {
  if (warnedTags.has(message))
    return;
  warnedTags.add(message);
  console.warn(`u2 spec: ${message}`);
}

/** The type rules the renderer enforces, as a query — null when `value` fits `prop`, the message
 * the renderer would throw otherwise. The designer validates an edit against the same code. */
export function checkProp(prop: SpecPropMeta, value: unknown): string | null {
  const type = prop.propertyType ?? prop.type;
  if (!fitsType(type, value))
    return `prop "${prop.name}" expects ${type}`;
  const choices = prop.choices;
  // scalar values only: on a string_list, choices constrain the items a picker offers, not the array
  if (choices && choices.length > 0 && type !== 'string_list' && type !== 'object' &&
    !choices.includes(value as string))
    return `prop "${prop.name}" must be one of ${choices.map((c) => `"${c}"`).join(', ')}`;
  return null;
}

function fitsType(type: string, value: unknown): boolean {
  switch (type) {
    case 'string':
      return typeof value === 'string';
    case 'int':
      return typeof value === 'number' && Number.isInteger(value);
    case 'double':
      return typeof value === 'number' && isFinite(value);
    case 'bool':
      return typeof value === 'boolean';
    case 'string_list':
      return Array.isArray(value) && value.every((item) => typeof item === 'string');
    case 'object':
      return isJson(value);
    default:
      return false;
  }
}

export function isHtmlTag(tag: string): boolean {
  return HTML_TAGS.has(tag);
}

/** The panel property model of a plain HTML node: text and cls everywhere, href on `a`,
 * src on `img` — exactly what the renderer honors. Everything is live-bindable (UB-4); the
 * appearance group rides by identity, which `sharedAppearance` and the manifest compare on. */
export function htmlProps(tag: string): SpecPropMeta[] {
  const props: SpecPropMeta[] = [
    {name: 'text', type: 'string', bindable: true, description: 'Text content'},
    {name: 'cls', type: 'string', bindable: true, description: 'Space-separated CSS classes'},
  ];
  if (tag === 'a')
    props.push({name: 'href', type: 'string', bindable: true, description: 'Link target'});
  if (tag === 'img')
    props.push({name: 'src', type: 'string', bindable: true, description: 'Image source'});
  props.push(...APPEARANCE_PROPS);
  return props;
}

function isJson(value: unknown): boolean {
  if (!JSON_TYPES.has(typeof value))
    return false;
  try {
    JSON.stringify(value, (_name, item) => {
      if (item === undefined || typeof item === 'function')
        throw new Error('not JSON');
      return item;
    });
    return true;
  } catch (e) {
    return false;
  }
}

/** Validates the envelope. Throws — a spec that is not a spec has no containable node. */
export function parseSpec(json: string | object): Spec {
  const parsed = typeof json === 'string' ? JSON.parse(json) : json;
  if (!parsed || typeof parsed !== 'object')
    throw new Error('u2 spec: expected a spec object or a JSON string');
  const spec = parsed as Partial<Spec>;
  if (spec.$schema !== SPEC_SCHEMA)
    throw new Error(`u2 spec: $schema must be "${SPEC_SCHEMA}", got ${JSON.stringify(spec.$schema)}`);
  if (!spec.root || typeof spec.root !== 'object' || typeof spec.root.tag !== 'string')
    throw new Error('u2 spec: "root" must be a node with a "tag"');
  if (spec.components !== undefined) {
    if (!Array.isArray(spec.components))
      throw new Error('u2 spec: "components" must be an array of nodes');
    for (const node of spec.components) {
      if (!node || typeof node !== 'object' || typeof node.tag !== 'string')
        throw new Error('u2 spec: every entry of "components" must be a node with a "tag"');
    }
  }
  return spec as Spec;
}

export function renderSpec(spec: string | object, ctx?: SpecContext, reg?: Registry,
  options?: SpecInstanceOptions): SpecInstance {
  return new SpecInstance(parseSpec(spec), ctx ?? new SpecContext(), reg, options);
}
