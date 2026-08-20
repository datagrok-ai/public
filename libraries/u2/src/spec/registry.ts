import {Control} from '../core/component.js';
import type {PropertyLike} from '../core/property-like.js';
import type {SpecNode} from './spec.js';

/** Envelope version of the JSON spec: `{"$schema": "dg-ui/1", "root": {…}}`. */
export const SPEC_SCHEMA = 'dg-ui/1';

/** A component prop, described the way the platform describes every property — so the property
 * grid, validation and the manifest all read one shape. `type` carries a platform TYPE string
 * ('string', 'int', 'double', 'bool', 'string_list', 'object'); `object` takes any
 * JSON-serializable payload — the escape hatch for structured options. */
export interface SpecPropMeta extends PropertyLike {
  type: string;
  /** Accepts a `bind` entry; the resolved context signal is passed to `create` as this prop. */
  bindable?: boolean;
  /** Edits made in the component flow back into the bound context signal. */
  twoWay?: boolean;
}

export interface ComponentMeta {
  /** `u2-*`. */
  tag: string;
  create: (props: Record<string, unknown>) => Control;
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
  /** Props a child node may carry that this component reads — an accordion pane's `title`. */
  childProps?: SpecPropMeta[];
  events?: string[];
  /** Spec children go into the component root; without it they are dropped with a warning. */
  acceptsChildren?: boolean;
  /** A valid node spec — the gallery renders it, an LLM learns the shape from it. */
  example: object;
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
    const withChildren = meta.createWithChildren;
    const stamped: ComponentMeta = {
      ...meta,
      create: (props) => Registry._stamp(meta.create(props), stamped, props),
      createWithChildren: withChildren ?
        (props, children, nodes) =>
          Registry._stamp(withChildren(props, children, nodes), stamped, props) :
        undefined,
    };
    this._metas.set(meta.tag, stamped);
  }

  get(tag: string): ComponentMeta | undefined {
    return this._metas.get(tag);
  }

  /** Everything a spec author, an LLM or a generated wrapper needs — the metadata minus its hooks. */
  manifest(): object {
    const components: Omit<ComponentMeta, 'create' | 'createWithChildren' | 'adopt'>[] = [];
    for (const meta of this._metas.values()) {
      const {create: _create, createWithChildren: _createWithChildren, adopt: _adopt, ...rest} = meta;
      components.push(rest);
    }
    return {schema: SPEC_SCHEMA, recipes: 'docs/recipes/', components};
  }

  private static _stamp(component: Control, meta: ComponentMeta,
    props: Record<string, unknown>): Control {
    component.meta = meta;
    component.specProps = props;
    return component;
  }
}

export const registry = new Registry();
