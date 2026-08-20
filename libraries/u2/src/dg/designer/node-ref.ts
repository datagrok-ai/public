/* What the designer selects and navigates: a handle on one node of a rendered spec, plus the two
   models derived from the node map — the structure tree and the breadcrumb path. Platform-free, so
   both are testable without the shell; only the handler and the view below need datagrok-api. */
import {Control} from '../../core/component.js';
import type {ComponentMetaLike} from '../../core/widget-like.js';
import type {TreeNode} from '../../components/tree.js';
import {registry as globalRegistry} from '../../spec/registry.js';
import type {SpecInstance, SpecNode} from '../../spec/spec.js';

const SEPARATOR = ' › ';
/** The class the spec renderer gives a node it could not build (`css/spec.css`). */
const ERROR_CLASS = 'u2-spec-error';

export interface SpecTree {
  roots: TreeNode<SpecNode>[];
  /** Dotted index per node — the tree's stable ids, and what {@link idPath} walks. */
  ids: Map<SpecNode, string>;
}

/** The designer's selection handle: what `grok.shell.o` carries, and what the node ObjectHandler
 * renders the context panel from. */
export class SpecNodeRef {
  constructor(readonly instance: SpecInstance, readonly node: SpecNode) {}

  /** What the spec built here — a component, a plain element, or the error placeholder. */
  built(): Control | HTMLElement | undefined {
    return this.instance.nodes().get(this.node);
  }

  element(): HTMLElement | undefined {
    const built = this.built();
    return built instanceof Control ? built.root : built;
  }

  /** The registry metadata behind the tag: the component's own stamp, or the registry itself for a
   * node that failed to build and has no component to ask. */
  meta(): ComponentMetaLike | undefined {
    const built = this.built();
    return (built instanceof Control ? built.meta : undefined) ?? globalRegistry.get(this.node.tag);
  }

  /** Why the node renders as a placeholder, or null when it built. */
  error(): string | null {
    const built = this.built();
    if (built === undefined || built instanceof Control || !built.classList.contains(ERROR_CLASS))
      return null;
    return built.textContent;
  }

  /** `root › details › nameInput` — the status bar's line and the panel's identity row. */
  path(): string {
    const parents = parentMap(this.instance);
    const parts: string[] = [];
    for (let at: SpecNode | undefined = this.node; at; at = parents.get(at))
      parts.unshift(nodeLabel(at));
    return parts.join(SEPARATOR);
  }
}

/** A node's spec identity where it has one, its tag otherwise. */
export function nodeLabel(node: SpecNode): string {
  return node.name ?? node.tag;
}

/** The structure tree over what a spec actually rendered — a child the renderer dropped (a node
 * under a component that takes none) is not in the map and not in the tree. */
export function specTree(instance: SpecInstance): SpecTree {
  const ids = new Map<SpecNode, string>();
  const build = (node: SpecNode, id: string): TreeNode<SpecNode> => {
    ids.set(node, id);
    const children = childrenOf(instance, node).map((child, i) => build(child, `${id}.${i}`));
    return {id, label: nodeLabel(node), data: node, children: children.length > 0 ? children : undefined};
  };
  return {roots: rootsOf(instance).map((node, i) => build(node, String(i))), ids};
}

/** Every id on the way to `id`, root first — `VirtualTree.expandPath`'s argument. */
export function idPath(id: string): string[] {
  const parts = id.split('.');
  return parts.map((_part, i) => parts.slice(0, i + 1).join('.'));
}

/** Nodes the renderer replaced with its error placeholder. */
export function brokenCount(instance: SpecInstance): number {
  let count = 0;
  for (const built of instance.nodes().values()) {
    if (!(built instanceof Control) && built.classList.contains(ERROR_CLASS))
      count++;
  }
  return count;
}

function childrenOf(instance: SpecInstance, node: SpecNode): SpecNode[] {
  const nodes = instance.nodes();
  return (node.children ?? []).filter((child) => nodes.has(child));
}

function rootsOf(instance: SpecInstance): SpecNode[] {
  const children = new Set<SpecNode>();
  for (const node of instance.nodes().keys()) {
    for (const child of childrenOf(instance, node))
      children.add(child);
  }
  const roots: SpecNode[] = [];
  for (const node of instance.nodes().keys()) {
    if (!children.has(node))
      roots.push(node);
  }
  return roots;
}

function parentMap(instance: SpecInstance): Map<SpecNode, SpecNode> {
  const parents = new Map<SpecNode, SpecNode>();
  for (const node of instance.nodes().keys()) {
    for (const child of childrenOf(instance, node))
      parents.set(child, node);
  }
  return parents;
}
