/* dg-ui/1: JSON in, a live component tree out, and back through dump(). A malformed node becomes
   a visible placeholder — the rest of the tree always renders. Specs never carry code: events
   name context commands, and that is the only executable form. */
import {signal, Signal} from '../core/signals.js';
import {Control} from '../core/component.js';
import {ComponentMeta, SpecPropMeta, Registry, registry as globalRegistry, SPEC_SCHEMA} from './registry.js';

export interface SpecNode {
  tag: string;
  /** Unique within the spec: what selection, patches, code-behind and automation address. */
  name?: string;
  props?: Record<string, unknown>;
  /** `{propName: '$.path'}` — binds a prop to a context signal. */
  bind?: Record<string, string>;
  /** `{eventName: 'cmd:name'}` — the only executable form a spec has. */
  on?: Record<string, string>;
  children?: SpecNode[];
}

export interface Spec {
  $schema: string;
  root: SpecNode;
}

export interface SpecContextOptions {
  data?: Record<string, unknown>;
  commands?: Record<string, () => void>;
}

const HTML_TAGS = new Set(['div', 'span', 'p', 'h1', 'h2', 'h3', 'a', 'img']);
const HTML_PROPS = new Set(['text', 'cls', 'href', 'src']);
const JSON_TYPES = new Set(['object', 'string', 'number', 'boolean']);

/** What `$.path` and `cmd:name` resolve against. Plain data values are wrapped into signals,
 * so a spec binds to them and sees live updates either way. */
export class SpecContext {
  readonly data: Record<string, Signal<unknown>> = {};
  readonly commands: Record<string, () => void>;

  constructor(options: SpecContextOptions = {}) {
    const data = options.data ?? {};
    for (const key of Object.keys(data)) {
      const value = data[key];
      this.data[key] = value instanceof Signal ? value as Signal<unknown> : signal(value);
    }
    this.commands = options.commands ?? {};
  }

  resolve(path: string): Signal<unknown> {
    const found = path.startsWith('$.') ? this.data[path.slice(2)] : undefined;
    if (!(found instanceof Signal))
      throw new Error(`nothing bound at "${path}"`);
    return found;
  }
}

/** A rendered spec: one scope owning every component the spec built, plus the node → built map
 * that makes {@link dump} read live state and {@link nodeAt} hit-test the DOM. */
export class SpecInstance extends Control {
  readonly ctx: SpecContext;

  private readonly _spec: Spec;
  private readonly _registry: Registry;
  private readonly _nodes = new Map<SpecNode, Control | HTMLElement>();
  private readonly _named = new Map<string, Control | HTMLElement>();
  private readonly _elements = new Map<Element, SpecNode>();
  private readonly _warned = new Set<string>();

  constructor(spec: Spec, ctx: SpecContext, reg: Registry = globalRegistry) {
    super();
    this._spec = spec;
    this.ctx = ctx;
    this._registry = reg;
    this.root.classList.add('u2-spec');
    this.run(() => this.root.append(SpecInstance._element(this._render(spec.root))));
  }

  dispose(): void {
    super.dispose();
    this._nodes.clear();
    this._named.clear();
    this._elements.clear();
  }

  /** The spec back out, with live values folded in: bound props keep their `bind` entry, plain
   * ones report what the component holds now. Rendering the result reproduces this UI. */
  dump(): object {
    return {$schema: SPEC_SCHEMA, root: this._dump(this._spec.root)};
  }

  /** Every node the spec rendered, components and plain HTML alike. */
  nodes(): ReadonlyMap<SpecNode, Control | HTMLElement> {
    return this._nodes;
  }

  node(name: string): Control | HTMLElement | undefined {
    return this._named.get(name);
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

  private _render(node: SpecNode, parent?: ComponentMeta): Control | HTMLElement {
    const built = this._build(node, parent);
    const el = SpecInstance._element(built);
    this._nodes.set(node, built);
    this._elements.set(el, node);
    if (node.name !== undefined) {
      // the automation id: what a test, a tutorial or an agent addresses the node by
      el.dataset.u2Name = node.name;
      if (this._named.has(node.name))
        this._warn(`${node.tag}: duplicate name "${node.name}"`);
      else
        this._named.set(node.name, built);
    }
    return built;
  }

  /** `parent` is the meta of the component this node hangs under: a prop the parent declares in
   * `childProps` — a pane title — validates against it instead of counting as unknown. */
  private _build(node: SpecNode, parent?: ComponentMeta): Control | HTMLElement {
    try {
      const commands = SpecInstance._commands(node);
      const meta = this._registry.get(node.tag);
      const built = meta ? this._component(node, meta, parent) : this._html(node, parent);
      for (const [event, command] of commands)
        this._listen(node, SpecInstance._element(built), event, command);
      return built;
    } catch (e) {
      return SpecInstance._placeholder(`${node.tag}: ${e instanceof Error ? e.message : String(e)}`);
    }
  }

  private _component(node: SpecNode, meta: ComponentMeta, parent?: ComponentMeta): Control {
    const props: Record<string, unknown> = {};
    for (const [name, value] of Object.entries(node.props ?? {})) {
      const prop = SpecInstance._prop(meta, name, parent);
      if (!prop)
        this._warn(`${node.tag}: unknown prop "${name}"`);
      props[name] = prop ? SpecInstance._checked(prop, value) : value;
    }

    const twoWay: [string, Signal<unknown>][] = [];
    for (const [name, path] of Object.entries(node.bind ?? {})) {
      const prop = SpecInstance._prop(meta, name);
      if (!prop || !prop.bindable)
        throw new Error(`prop "${name}" is not bindable`);
      const source = this.ctx.resolve(path);
      props[name] = source;
      if (prop.twoWay)
        twoWay.push([name, source]);
    }

    const component = meta.createWithChildren ?
      meta.createWithChildren(props, this._children(node, meta, true), node.children ?? []) :
      this._adopt(node, meta, meta.create(props));
    for (const [name, source] of twoWay)
      this._bridge(component, name, source);
    return component;
  }

  private _html(node: SpecNode, parent?: ComponentMeta): HTMLElement {
    if (!HTML_TAGS.has(node.tag))
      throw new Error('neither a registered component nor a supported HTML tag');
    if (node.bind)
      throw new Error('only registered components can be bound');
    const el = document.createElement(node.tag);
    for (const [name, value] of Object.entries(node.props ?? {})) {
      if (!HTML_PROPS.has(name)) {
        if (!SpecInstance._childProp(parent, name))
          this._warn(`${node.tag}: unknown prop "${name}"`);
        continue;
      }
      if (typeof value !== 'string')
        throw new Error(`prop "${name}" expects string`);
      if (name === 'text')
        el.textContent = value;
      else if (name === 'cls')
        el.classList.add(...value.split(' ').filter((c) => c));
      else
        el.setAttribute(name, value);
    }
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

  /** Both directions, echo-suppressed by comparison: a write made inside an effect is flushed
   * after that effect returns, so a transient flag would already be down (property-grid's lesson).
   * Components that took the context signal itself need no bridge at all. */
  private _bridge(component: Control, name: string, source: Signal<unknown>): void {
    const own = (component as unknown as Record<string, unknown>)[name];
    if (!(own instanceof Signal) || own === source)
      return;
    this.scope.effect(() => {
      const value = source.value;
      if (own.peek() !== value)
        own.value = value;
    });
    this.scope.effect(() => {
      const value = own.value;
      if (source.peek() !== value)
        source.value = value;
    });
  }

  private _listen(node: SpecNode, el: HTMLElement, event: string, command: string): void {
    const handler = () => {
      const run = this.ctx.commands[command];
      if (typeof run !== 'function')
        this._warn(`${node.tag}: no command "${command}"`);
      else
        run();
    };
    el.addEventListener(event, handler);
    this.scope.own(() => el.removeEventListener(event, handler));
  }

  private _dump(node: SpecNode): object {
    const out: Record<string, unknown> = {tag: node.tag};
    if (node.name !== undefined)
      out.name = node.name;
    if (node.props)
      out.props = this._dumpProps(node);
    if (node.bind)
      out.bind = {...node.bind};
    if (node.on)
      out.on = {...node.on};
    if (node.children)
      out.children = node.children.map((child) => this._dump(child));
    return out;
  }

  private _dumpProps(node: SpecNode): Record<string, unknown> {
    const live = this._liveValue(node);
    const props: Record<string, unknown> = {};
    for (const [name, value] of Object.entries(node.props ?? {}))
      props[name] = live && name === 'value' ? live.current : value;
    return props;
  }

  /** The component's own `value` signal, but only where the spec's prop — not a binding — feeds it. */
  private _liveValue(node: SpecNode): {current: unknown} | null {
    const meta = this._registry.get(node.tag);
    const component = this._nodes.get(node);
    if (!meta || !(component instanceof Control) || (node.bind && 'value' in node.bind))
      return null;
    if (!meta.props.some((p) => p.name === 'value'))
      return null;
    const own = (component as unknown as {value?: unknown}).value;
    return own instanceof Signal ? {current: own.peek()} : null;
  }

  private _warn(message: string): void {
    if (this._warned.has(message))
      return;
    this._warned.add(message);
    console.warn(`u2 spec: ${message}`);
  }

  private static _commands(node: SpecNode): [string, string][] {
    const commands: [string, string][] = [];
    for (const [event, command] of Object.entries(node.on ?? {})) {
      if (!command.startsWith('cmd:'))
        throw new Error(`on.${event} must be "cmd:<name>" — a spec never carries code`);
      commands.push([event, command.slice(4)]);
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
    return built instanceof Control ? built.root : built;
  }

  private static _checked(prop: SpecPropMeta, value: unknown): unknown {
    const type = prop.propertyType ?? prop.type;
    switch (type) {
      case 'string':
        if (typeof value === 'string')
          return value;
        break;
      case 'int':
        if (typeof value === 'number' && Number.isInteger(value))
          return value;
        break;
      case 'double':
        if (typeof value === 'number' && isFinite(value))
          return value;
        break;
      case 'bool':
        if (typeof value === 'boolean')
          return value;
        break;
      case 'string_list':
        if (Array.isArray(value) && value.every((item) => typeof item === 'string'))
          return value;
        break;
      case 'object':
        if (SpecInstance._isJson(value))
          return value;
        break;
    }
    throw new Error(`prop "${prop.name}" expects ${type}`);
  }

  private static _isJson(value: unknown): boolean {
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

  private static _placeholder(message: string): HTMLElement {
    const el = document.createElement('div');
    el.className = 'u2-spec-error';
    el.textContent = message;
    return el;
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
  return spec as Spec;
}

export function renderSpec(spec: string | object, ctx?: SpecContext, reg?: Registry): SpecInstance {
  return new SpecInstance(parseSpec(spec), ctx ?? new SpecContext(), reg);
}
