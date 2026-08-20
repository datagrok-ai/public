/* The binding picker's model (Q7): everything a node may bind to — context data, the tray, and the
   named controls — as one tree walked through the DD5 BindSource protocol. Pure and testable
   without the platform: the picker renders it, the headless suite reads it. Enumeration allocates
   nothing (the laziness contract), so every branch is a thunk read when the user opens it. */
import {Component} from '../../core/component.js';
import {isBindSource} from '../../spec/bind-source.js';
import type {BindProp, BindSource} from '../../spec/bind-source.js';
import {segmentsToPath} from '../../spec/path.js';
import type {SpecInstance, SpecNode} from '../../spec/spec.js';

export interface BindTreeNode {
  label: string;
  /** What picking this node binds to; null for a group — a walkable step that is nothing on its
   * own (`currentRow`). */
  path: string | null;
  prop?: BindProp;
  children?: () => BindTreeNode[];
}

/** The roots, in the order the picker offers them: context data, the tray, then the named controls
 * that carry signal-backed props. A source's own steps hang under it. */
export function bindTree(instance: SpecInstance): BindTreeNode[] {
  const roots: BindTreeNode[] = [];
  for (const prop of instance.ctx.bindProps()) {
    const value = instance.ctx.data[prop.name];
    roots.push(root(labelOf(prop), prop.name, isBindSource(value) ? value : undefined, prop));
  }
  for (const node of instance.spec.components ?? []) {
    const source = sourceOf(instance, node);
    if (source !== undefined)
      roots.push(root(`${node.name} (${node.tag})`, node.name!, source));
  }
  for (const node of named(instance.spec.root)) {
    const source = sourceOf(instance, node);
    if (source !== undefined && source.bindProps().length > 0)
      roots.push(root(`${node.name} (${node.tag})`, node.name!, source));
  }
  return roots;
}

/** A broken component maps the placeholder ELEMENT, not a component — the guard is what keeps the
 * tray's failures out of the picker instead of throwing at it. */
function sourceOf(instance: SpecInstance, node: SpecNode): BindSource | undefined {
  const built = instance.nodes().get(node);
  return node.name !== undefined && built instanceof Component && isBindSource(built) ? built : undefined;
}

/** What an expander says when it has nothing to offer YET — a source that has not run knows no
 * columns, and an empty branch is indistinguishable from a broken one. Not pickable (no path). */
export const NOTHING_YET = 'nothing here yet — Refresh the source, or set designData to sample';

function root(label: string, name: string, source: BindSource | undefined,
  prop?: BindProp): BindTreeNode {
  const path = segmentsToPath([name]);
  return source === undefined ? {label, path, prop} :
    {label, path, prop, children: () => stepsOf(source, [name])};
}

function stepsOf(source: BindSource, segments: string[]): BindTreeNode[] {
  const steps = source.bindProps().map((prop): BindTreeNode => {
    const at = [...segments, prop.name];
    if (prop.walkable !== true)
      return {label: labelOf(prop), path: segmentsToPath(at), prop};
    return {label: labelOf(prop), path: null, prop, children: () => {
      const next = source.bindStep(prop.name);
      return isBindSource(next) ? stepsOf(next, at) : [];
    }};
  });
  return steps.length > 0 ? steps : [{label: NOTHING_YET, path: null}];
}

/** Every named node of the form, in document order. */
function named(node: SpecNode): SpecNode[] {
  const found = node.name === undefined ? [] : [node];
  for (const child of node.children ?? [])
    found.push(...named(child));
  return found;
}

function labelOf(prop: BindProp): string {
  const type = prop.propertyType ?? prop.type;
  let label = prop.name;
  if (type != null && type !== '')
    label += ` : ${type}`;
  if (prop.semType != null && prop.semType !== '')
    label += ` · ${prop.semType}`;
  // the two-way affordance: this leaf takes edits back
  return prop.writable === true ? `${label} ⇄` : label;
}
