import {Component, Control} from '../core/component.js';
import type {Signal} from '../core/signals.js';
import {APPEARANCE_PROPS, appearanceFor, applyAppearance} from './appearance.js';
import type {NamedProperty} from '../core/widget-like.js';
import type {ComponentMetaBase} from '../core/widget-like.js';
import type {SpecInstance, SpecNode} from './spec.js';

/** Envelope version of the JSON spec: `{"$schema": "dg-ui/1", "root": {…}}`. */
export const SPEC_SCHEMA = 'dg-ui/1';

/** The binding rule, stated once at the manifest's top level (UB-9) — no per-prop repetition. */
const BINDING_RULE = 'Every declared prop accepts a "bind" entry and follows it live, except props ' +
  'marked structural, which read a bind once. A bind on an undeclared prop is an error.';

/** Injected wherever the appearance group is (every visual tag but a viewer, whose look is its own),
 * listed with the component's own props; live through `Control.visible`. */
const VISIBLE: SpecPropMeta = {name: 'visible', type: 'bool', description: 'false hides the node; true by default.'};

/** A component prop, described the way the platform describes every property — so the property
 * grid, validation and the manifest all read one shape. `type` carries a platform TYPE string
 * ('string', 'int', 'double', 'bool', 'string_list', 'object'); `object` takes any
 * JSON-serializable payload — the escape hatch for structured options. */
export interface SpecPropMeta extends NamedProperty {
  type: string;
  /** Accepts a `bind` entry; the resolved context signal is passed to `create` as this prop. */
  bindable?: boolean;
  /** Edits made in the component flow back into the bound context signal. */
  twoWay?: boolean;
  /** An object prop that accepts dotted `bind` keys (`params.days`) — sub-binds, delivered raw
   * through the components engine and resolved at start. */
  subBindable?: boolean;
  /** Only a rebuild can change it (an input's `inline`, a splitter's `direction`): a bind reads it
   * once and warns; the designer's set-prop patch rebuilds the node. Every other prop is live —
   * bindable, or written through the component's setter on every change. */
  structural?: boolean;
}

/** What a design-time action wants done, described against the node's current state — the
 * designer view converts it into an engine patch. */
export type DesignerActionPatch =
  | {op: 'add-child', node: SpecNode}
  | {op: 'set-prop', name: string, value: unknown};

/** A per-component design-time verb ('Add tab'), shown on the designer's context menus. */
export interface DesignerAction {
  /** The verb's caption in the designer's menus. */
  text: string;
  /** A platform icon name ('filter') for the menus that show it. */
  icon?: string;
  /** Reads the node's current children/props — and what it built, when rendered — and answers
   * what to do; null is nothing to do (cancelled). May ask the user first, hence the promise. */
  produce: (node: SpecNode, built?: Component | HTMLElement) =>
    DesignerActionPatch | null | Promise<DesignerActionPatch | null>;
}

/** What a non-visual component is built with (Q9): the mode it builds for, and the binds it
 * resolves once the whole tree exists. */
export interface ComponentEnv {
  designTime: boolean;
  /** Resolution deferred to {@link ComponentStart.start}: by then a param may bind to an input
   * declared after this source. Null where the path does not resolve — the renderer warns. */
  resolve: (path: string) => Signal<unknown> | null;
  /** The node's raw `bind` entries, dotted sub-bind keys included ('params.days' → '$.daysInput'). */
  subBinds: Record<string, string>;
}

/** The second construction phase (Q1) a component may declare: run once the visual tree exists,
 * where an auto-run kicks and a sub-bind resolves. */
export interface ComponentStart {
  start(): void;
}

export interface ComponentMeta extends ComponentMetaBase {
  /** `u2-*`. */
  tag: string;
  /** What the palette shows for the tag ('Scatter plot'); the tag's suffix where absent. */
  label?: string;
  /** Palette and manifest grouping: 'Inputs' | 'Containers' | 'Actions' | 'Display' | 'Data'. */
  category?: string;
  /** False puts the tag on the non-visual tray: a data source, built through
   * {@link createComponent} and never appended to the canvas. */
  visual?: boolean;
  /** False keeps the shared Appearance prop block off this tag — viewers, whose look has its own
   * color semantics. */
  appearance?: boolean;
  create?: (props: Record<string, unknown>) => Control;
  /** How a `visual: false` tag is built — no DOM, and the env carries the design-time mode and
   * the binds it resolves at {@link ComponentStart.start}. */
  createComponent?: (props: Record<string, unknown>, env: ComponentEnv) => Component;
  /** Replaces {@link create} where children are constructor arguments (splitter panels, tabs).
   * They arrive rendered, with their spec nodes alongside so per-child props can be read. */
  createWithChildren?: (props: Record<string, unknown>, children: (Control | HTMLElement)[],
    nodes: SpecNode[]) => Control;
  /** Hands each rendered child to the component's own add-method (`Form.add`) instead of appending
   * it to the root; implies children are accepted. */
  adopt?: (parent: Control, child: Control | HTMLElement, index: number) => void;
  /** Moves view state a designer rebuild would lose (a tab strip's active tab) from the outgoing
   * instance to the fresh one; a closure, so {@link Registry.manifest} strips it. */
  carry?: (node: SpecNode, prev: Component | HTMLElement, next: Component | HTMLElement) => void;
  /** Writes a bound {@link childProps} entry of the child at `index` in place (a tab's title) on
   * every change of its source — every child prop not marked structural. A closure, so
   * {@link Registry.manifest} strips it. */
  setChildProp?: (parent: Control, index: number, name: string, value: unknown) => void;
  description: string;
  /** A sentence of judgment for the manifest reader (an LLM or a reviewer): when to reach for
   * this component, and when something else serves better. Reference facts stay in
   * {@link description}; composition-level judgment lives in docs/recipes. */
  usage?: string;
  props: SpecPropMeta[];
  /** Props a designer seeds a freshly inserted node with ('Item 1..3' for breadcrumbs) —
   * plain JSON, so it stays in the manifest. */
  defaults?: Record<string, unknown>;
  /** Children a designer seeds a freshly inserted node with — a multi-host (tabs, accordion) arrives
   * with one pane, so there is something to see and click. Deep-cloned per insert and named like
   * any inserted node; plain JSON, so it stays in the manifest. */
  defaultChildren?: SpecNode[];
  /** Props a child node may carry that this component reads — an accordion pane's `title`. */
  childProps?: SpecPropMeta[];
  events?: string[];
  /** Spec children go into the component root; without it they are dropped with a warning. */
  acceptsChildren?: boolean;
  /** A valid node spec — the gallery renders it, an LLM learns the shape from it. */
  example: object;
  /** The tag-level sample a design-time preview falls back to where the node carries no `sample`
   * prop (DD9); plain JSON, so it stays in the manifest. */
  designPreview?: object;
  /** Carries closures, so {@link Registry.manifest} strips it. */
  designerActions?: DesignerAction[];
  /** A palette glyph (the platform's own viewer icon); a closure, so {@link Registry.manifest} strips it. */
  icon?: () => HTMLElement;
}

/** A companion library (ems, later chem/bio control sets) plugging into a registry through
 * {@link Registry.use}. */
export interface U2Plugin {
  /** One install per name and registry: a second copy of the library on the page joins the first. */
  name: string;
  /** Tag prefixes the plugin owns ('u2-domain-'): an unregistered tag with one names the plugin. */
  tagPrefixes?: string[];
  register(registry: Registry): void;
  /** Wraps the build of every tray component of a spec rendered over the registry — the ems unit
   * of work. Several plugins nest in install order, the first installed outermost. */
  specAmbient?: <T>(instance: SpecInstance, build: () => T) => T;
}

/** Libraries whose tags a stored spec may name before the host has installed them. */
const KNOWN_PLUGINS: Record<string, string> = {
  'u2-domain-': '@datagrok-libraries/ems — call installEms()',
  'u2-manifest-': '@datagrok-libraries/ems — call installEms()',
};

/** Tag → component metadata, and the source of the machine-readable manifest. */
export class Registry {
  private readonly _metas = new Map<string, ComponentMeta>();
  private readonly _plugins: {plugin: U2Plugin, uninstall: () => void}[] = [];
  /** Set by `registerAll()`: u2's own components are registered here. */
  builtins = false;

  /** Installs a plugin once per name: its tags register, and its spec ambient wraps every spec
   * rendered over this registry from then on. Answers the uninstall, which removes the tags it
   * registered and the ambient; a second install of the same name answers the first one's. */
  use(plugin: U2Plugin): () => void {
    const installed = this._plugins.find((p) => p.plugin.name === plugin.name);
    if (installed !== undefined) {
      if (installed.plugin !== plugin)
        console.warn(`u2 registry: "${plugin.name}" is installed already, from another copy — that one stays`);
      return installed.uninstall;
    }
    const before = new Set(this._metas.keys());
    plugin.register(this);
    const tags = [...this._metas.keys()].filter((tag) => !before.has(tag));
    const entry = {plugin, uninstall: (): void => {
      const at = this._plugins.indexOf(entry);
      if (at < 0)
        return;
      this._plugins.splice(at, 1);
      for (const tag of tags)
        this._metas.delete(tag);
    }};
    this._plugins.push(entry);
    return entry.uninstall;
  }

  /** Install order — the order the spec ambients nest in. */
  plugins(): readonly U2Plugin[] {
    return this._plugins.map((p) => p.plugin);
  }

  /** What would register an unknown tag: the installed plugin owning its prefix, else the library
   * known to own it; null when nothing claims the tag. */
  missingTagHint(tag: string): string | null {
    const owner = this._plugins.find((p) => p.plugin.tagPrefixes?.some((prefix) => tag.startsWith(prefix)));
    if (owner !== undefined)
      return `${owner.plugin.name} is installed and registers no such tag`;
    const known = Object.keys(KNOWN_PLUGINS).find((prefix) => tag.startsWith(prefix));
    if (known !== undefined)
      return `it comes with ${KNOWN_PLUGINS[known]}`;
    // u2's own components are registered by a call, not by importing the library
    return tag.startsWith('u2-') && !this.builtins ?
      'u2\'s components are not registered in this registry — call registerAll() first' : null;
  }

  /** The construction hooks are wrapped so that everything the registry builds carries its
   * metadata and the props it was built from — the source of the component's introspection
   * surface. A component constructed by hand (`new TextInput()`) carries neither, and answers an
   * empty property list. */
  register(meta: ComponentMeta): void {
    if (this._metas.has(meta.tag))
      throw new Error(`u2 registry: "${meta.tag}" is already registered`);
    const create = meta.create;
    const withChildren = meta.createWithChildren;
    const component = meta.createComponent;
    if (meta.visual === false ? component === undefined : create === undefined && withChildren === undefined) {
      throw new Error(`u2 registry: "${meta.tag}" needs ${meta.visual === false ?
        'createComponent' : 'create or createWithChildren'}`);
    }
    const looks = meta.visual !== false && meta.appearance !== false;
    const injected = looks ? appearanceFor(meta.props) : [];
    const visible = looks && !meta.props.some((p) => p.name === VISIBLE.name) ? [VISIBLE] : [];
    const finish = (built: Control, props: Record<string, unknown>): Control => {
      Registry._stamp(built, stamped, props);
      applyAppearance(built, props, injected);
      if (visible.length > 0 && typeof props.visible === 'boolean')
        built.visible = props.visible;
      return built;
    };
    const stamped: ComponentMeta = {
      ...meta,
      props: [...meta.props, ...visible, ...injected],
      create: create ? (props) => finish(create(props), props) : undefined,
      createWithChildren: withChildren ?
        (props, children, nodes) => finish(withChildren(props, children, nodes), props) :
        undefined,
      createComponent: component ?
        (props, env) => Registry._stamp(component(props, env), stamped, props) :
        undefined,
    };
    this._metas.set(meta.tag, stamped);
  }

  get(tag: string): ComponentMeta | undefined {
    return this._metas.get(tag);
  }

  /** Registration order — the palette's source. */
  metas(): readonly ComponentMeta[] {
    return [...this._metas.values()];
  }

  /** Everything a spec author, an LLM or a generated wrapper needs — the metadata minus its hooks. */
  manifest(): object {
    const components: Omit<ComponentMeta, 'create' | 'createWithChildren' | 'createComponent' | 'adopt' |
      'carry' | 'setChildProp' | 'designerActions' | 'icon'>[] = [];
    for (const meta of this._metas.values()) {
      const {create: _create, createWithChildren: _createWithChildren, createComponent: _createComponent,
        adopt: _adopt, carry: _carry, setChildProp: _setChildProp, designerActions: _designerActions,
        icon: _icon, ...rest} = meta;
      rest.props = rest.props.filter((p) => !APPEARANCE_PROPS.includes(p));
      components.push(rest);
    }
    return {schema: SPEC_SCHEMA, recipes: 'docs/recipes/', binding: BINDING_RULE,
      appearance: APPEARANCE_PROPS, components};
  }

  private static _stamp<T extends Component>(component: T, meta: ComponentMeta,
    props: Record<string, unknown>): T {
    component.componentMeta = meta;
    component.specProps = props;
    return component;
  }
}

export const registry = new Registry();
