import {Component} from '../core/component.js';
import type {SpecNode} from './spec.js';

/** Envelope version of the JSON spec: `{"$schema": "dg-ui/1", "root": {…}}`. */
export const SPEC_SCHEMA = 'dg-ui/1';

export interface PropMeta {
  name: string;
  /** `json` takes any JSON-serializable payload — the escape hatch for structured options. */
  type: 'string' | 'int' | 'float' | 'bool' | 'string[]' | 'json';
  /** Accepts a `bind` entry; the resolved context signal is passed to `create` as this prop. */
  bindable?: boolean;
  /** Edits made in the component flow back into the bound context signal. */
  twoWay?: boolean;
  description?: string;
}

export interface ComponentMeta {
  /** `u2-*`. */
  tag: string;
  create: (props: Record<string, unknown>) => Component;
  /** Replaces {@link create} where children are constructor arguments (splitter panels, tabs).
   * They arrive rendered, with their spec nodes alongside so per-child props can be read. */
  createWithChildren?: (props: Record<string, unknown>, children: (Component | HTMLElement)[],
    nodes: SpecNode[]) => Component;
  /** Hands each rendered child to the component's own add-method (`Form.add`) instead of appending
   * it to the root; implies children are accepted. */
  adopt?: (parent: Component, child: Component | HTMLElement, index: number) => void;
  description: string;
  /** A sentence of judgment for the manifest reader (an LLM or a reviewer): when to reach for
   * this component, and when something else serves better. Reference facts stay in
   * {@link description}; composition-level judgment lives in docs/recipes. */
  usage?: string;
  props: PropMeta[];
  /** Props a child node may carry that this component reads — an accordion pane's `title`. */
  childProps?: PropMeta[];
  events?: string[];
  /** Spec children go into the component root; without it they are dropped with a warning. */
  acceptsChildren?: boolean;
  /** A valid node spec — the gallery renders it, an LLM learns the shape from it. */
  example: object;
}

/** Tag → component metadata, and the source of the machine-readable manifest. */
export class Registry {
  private readonly _metas = new Map<string, ComponentMeta>();

  register(meta: ComponentMeta): void {
    if (this._metas.has(meta.tag))
      throw new Error(`u2 registry: "${meta.tag}" is already registered`);
    this._metas.set(meta.tag, meta);
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
}

export const registry = new Registry();
