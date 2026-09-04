/* The binding picker (Q7): the bind tree in a dialog, grouped by where each root came from, with one
   search box over the labels. It answers the assembled path and nothing else — the commit stays with
   the panel, whose funnel (canApply → set-bind → warn) is the one authority on what may be written. */
import * as grok from 'datagrok-api/grok';
import {Dialog} from '../../components/containers/dialog.js';
import {VirtualTree} from '../../components/collections/tree.js';
import type {TreeNode} from '../../components/collections/tree.js';
import {TextInput} from '../../components/inputs/text-input.js';
import {iconButton} from '../../components/actions/buttons.js';
import {div, divV, span} from '../../core/elements.js';
import type {Input} from '../../core/input-base.js';
import {parsePath} from '../../spec/path.js';
import type {SpecInstance} from '../../spec/spec.js';
import {bindTree, named} from './bind-model.js';
import type {BindTreeNode} from './bind-model.js';

/** One heading of the picker: the roots under it, and where they came from — a form that offers
   `reagent` before the user has made anything has to say who put it there. */
export interface BindGroup {
  title: string;
  nodes: BindTreeNode[];
}

/** What a query shows, and the branches that have to be open for its matches to be visible. */
export interface BindRows {
  rows: TreeNode<BindTreeNode>[];
  expand: string[];
}

const APP_DATA = 'App data';
const SOURCES = 'Data sources';
const VIEWERS = 'Viewers';
const CONTROLS = 'Form controls';
/** The picker's one rule, as its hint and as the warning OK gives when a group is what is selected. */
export const PICK_A_VALUE = 'Pick a value — a group has nothing to bind to on its own';

/** The roots under the heading that says where they came from. A component wins a name collision
 * with a context key (OR-2): the spec's own declaration is the closer one. A platform viewer
 * (registry category `Viewers`) is a frame to walk, not a control to follow — its own heading. */
export function bindGroups(instance: SpecInstance, roots: BindTreeNode[]): BindGroup[] {
  const sources = new Set<string>();
  for (const node of instance.spec.components ?? []) {
    if (node.name !== undefined)
      sources.add(node.name);
  }
  const viewers = new Set<string>();
  for (const node of named(instance.spec.root)) {
    if (instance.registry.get(node.tag)?.category === VIEWERS)
      viewers.add(node.name!);
  }
  const data = new Set(Object.keys(instance.ctx.data));
  const groups = new Map<string, BindTreeNode[]>([[APP_DATA, []], [SOURCES, []], [VIEWERS, []], [CONTROLS, []]]);
  for (const root of roots) {
    const name = root.name ?? '';
    groups.get(sources.has(name) ? SOURCES : data.has(name) ? APP_DATA : viewers.has(name) ? VIEWERS : CONTROLS)!
      .push(root);
  }
  return [...groups].filter(([, nodes]) => nodes.length > 0).map(([title, nodes]) => ({title, nodes}));
}

/** The rows for `query`, and the branches that reveal its matches. An id is the node's own place in
 * the tree rather than the query it was built under, so an expansion survives the next keystroke —
 * the search used to reset every one of them. A match keeps its lazy branch, so its whole subtree
 * still browses; a branch kept only because a descendant matched carries its filtered children as an
 * array, which the tree never caches under the id the unfiltered branch will want back. */
export function bindRows(nodes: BindTreeNode[], query: string, prefix = '', depth = 3): BindRows {
  const text = query.trim().toLowerCase();
  const rows: TreeNode<BindTreeNode>[] = [];
  const expand: string[] = [];
  for (const node of nodes) {
    const id = `${prefix}${node.label}`;
    const children = node.children;
    if (text === '' || node.label.toLowerCase().includes(text)) {
      rows.push({id, label: node.label, data: node, children: children === undefined ? undefined :
        () => Promise.resolve(bindRows(children(), '', `${id}/`).rows)});
      continue;
    }
    const below = depth > 0 && children !== undefined ?
      bindRows(children(), text, `${id}/`, depth - 1) : undefined;
    if (below === undefined || below.rows.length === 0)
      continue;
    expand.push(id, ...below.expand);
    rows.push({id, label: node.label, data: node, children: below.rows});
  }
  return {rows, expand};
}

/** Who the picker is binding: the property, for the dialog's title, and the path it follows
 * today, whose row a re-pick arrives with revealed and selected. */
export interface BindPickerOptions {
  prop?: string;
  current?: string;
}

/** Opens the picker over everything `instance` can bind to. `onPick` gets the assembled path when
 * OK closes it on a leaf — answering `false` (a refused commit, balloon already up) keeps the
 * dialog open for a correction; on a group — a step that is nothing on its own — OK warns and the
 * dialog stays open, because a close that committed nothing read as a pick that silently failed.
 * The dialog is adopted by the ambient scope, so the panel that opened it takes it down with it. */
export function bindPicker(instance: SpecInstance, onPick: (path: string) => boolean | void,
  options: BindPickerOptions = {}): Dialog {
  const dialog = new Dialog(options.prop === undefined ? 'Bind to' : `${options.prop} — bind to`);
  dialog.root.classList.add('u2-bind-picker');
  const groups = bindGroups(instance, bindTree(instance));
  const search = dialog.runInScope(() =>
    new TextInput({search: true, inline: true, placeholder: 'Search bindings'}));
  const tree = dialog.runInScope(() => new VirtualTree<BindTreeNode>());
  const empty = span('', 'u2-picker-empty');
  tree.expanded.value = new Set(groups.map((group) => group.title));
  dialog.effect(() => {
    const query = search.value.value;
    const rows: TreeNode<BindTreeNode>[] = [];
    const expand: string[] = [];
    for (const group of groups) {
      const found = bindRows(group.nodes, query, `${group.title}/`);
      if (found.rows.length === 0)
        continue;
      rows.push({id: group.title, label: group.title, data: {label: group.title, path: null},
        children: found.rows});
      expand.push(...found.expand);
    }
    tree.setRoots(rows);
    // what the query matched is revealed, not merely kept: the match used to arrive collapsed
    // under three closed branches, which is what made the search worse than scrolling
    if (expand.length > 0)
      tree.expanded.value = new Set([...tree.expanded.peek(), ...expand]);
    const found = rows.length > 0;
    empty.textContent = found ? '' : query.trim() === '' ?
      'Nothing to bind to yet — add a data source, or name a control.' :
      `No bindings match "${query.trim()}".`;
    tree.root.style.display = found ? '' : 'none';
    empty.style.display = found ? 'none' : '';
  });

  dialog.add(divV([search.root, tree.root, empty,
    span(`${PICK_A_VALUE}. App data comes from the app that opened the designer.`, 'u2-bind-picker-hint')],
  'u2-bind-picker-body'));
  if (options.current !== undefined)
    void tree.expandPath(idsOf(groups, options.current));
  dialog.onCancel(() => dialog.dispose());
  dialog.onOK(() => {
    const picked = tree.selectedNode.peek()?.data;
    if (picked?.path == null) {
      grok.shell.warning(PICK_A_VALUE);
      return false;
    }
    if (onPick(picked.path) === false)
      return false;
    dialog.dispose();
  });
  return dialog.show({modal: true, width: 460});
}

/** The tree ids from a group heading down to `path`'s row — what a re-pick reveals. A step the
 * tree does not know (a source that has not run) ends the chain there, revealed as far as it goes.
 * The id scheme (`group/label/label…`) is assembled both here and in {@link bindRows} — keep the
 * two in sync. */
function idsOf(groups: BindGroup[], path: string): string[] {
  const segments = parsePath(path);
  if (segments === null)
    return [];
  const group = groups.find((g) => g.nodes.some((root) => root.name === segments[0]));
  if (group === undefined)
    return [];
  let node = group.nodes.find((root) => root.name === segments[0])!;
  const ids = [group.title, `${group.title}/${node.label}`];
  for (const segment of segments.slice(1)) {
    const next = node.children?.().find((child) => child.prop?.name === segment);
    if (next === undefined)
      break;
    ids.push(`${ids[ids.length - 1]}/${next.label}`);
    node = next;
  }
  return ids;
}

/** The `…` affordance beside a path field, on the input's options rail. `Input.addOptions` is
 * protected, so the rail is materialized the way the platform pickers materialize it. */
export function bindPickerButton(input: Input<any>, name: string, open: () => void,
  tooltip = 'Pick a binding'): void {
  let rail = input.box.querySelector('.u2-input-options') as HTMLElement | null;
  if (rail === null) {
    rail = div([], 'u2-input-options');
    input.box.append(rail);
  }
  const button = iconButton('ellipsis-h', open, {tooltip});
  button.dataset.u2BindPick = name;
  rail.append(button);
}
