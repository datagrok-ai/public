import {Component, Control} from '../core/component.js';
import type {Signal} from '../core/signals.js';
import {APPEARANCE_PROPS, appearanceFor, applyAppearance} from './appearance.js';
import type {NamedProperty} from '../core/widget-like.js';
import type {ComponentMetaBase} from '../core/widget-like.js';
import type {SpecNode} from './spec.js';

/** Envelope version of the JSON spec: `{"$schema": "dg-ui/1", "root": {…}}`. */
export const SPEC_SCHEMA = 'dg-ui/1';

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
}

/** What a design-time action wants done, described against the node's current state — the
 * designer view converts it into an engine patch. */
export type DesignerActionPatch =
  | {op: 'add-child', node: SpecNode}
  | {op: 'set-prop', name: string, value: unknown};

/** A per-component design-time verb ('Add tab'), shown on the designer's context menus. */
export interface DesignerAction {
  name: string;
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

/** Tag → component metadata, and the source of the machine-readable manifest. */
export class Registry {
  private readonly _metas = new Map<string, ComponentMeta>();

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
    const injected = meta.visual === false || meta.appearance === false ? [] : appearanceFor(meta.props);
    const stamped: ComponentMeta = {
      ...meta,
      props: [...meta.props, ...injected],
      create: create ? (props) => {
        const built = Registry._stamp(create(props), stamped, props);
        applyAppearance(built, props, injected);
        return built;
      } : undefined,
      createWithChildren: withChildren ?
        (props, children, nodes) => {
          const built = Registry._stamp(withChildren(props, children, nodes), stamped, props);
          applyAppearance(built, props, injected);
          return built;
        } :
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
    const components: Omit<ComponentMeta,
      'create' | 'createWithChildren' | 'createComponent' | 'adopt' | 'designerActions' | 'icon'>[] = [];
    for (const meta of this._metas.values()) {
      const {create: _create, createWithChildren: _createWithChildren, createComponent: _createComponent,
        adopt: _adopt, designerActions: _designerActions, icon: _icon, ...rest} = meta;
      rest.props = rest.props.filter((p) => !APPEARANCE_PROPS.includes(p));
      components.push(rest);
    }
    return {schema: SPEC_SCHEMA, recipes: 'docs/recipes/', appearance: APPEARANCE_PROPS, components};
  }

  private static _stamp<T extends Component>(component: T, meta: ComponentMeta,
    props: Record<string, unknown>): T {
    component.componentMeta = meta;
    component.specProps = props;
    return component;
  }
}

export const registry = new Registry();
