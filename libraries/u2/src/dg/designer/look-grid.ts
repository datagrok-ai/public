/* The viewer node's Properties section (VP-21..VP-24): the platform's own property grid over the
   live look — the viewer's categories, labels, column pickers — with every design-time commit
   captured from the viewer's property event as a `set-prop` patch. The grid follows the node's
   build: a patch re-creates the viewer, and the grid is re-pointed at the new look. */
import * as DG from 'datagrok-api/dg';
import {Component} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {span} from '../../core/elements.js';
import type {SpecPatch} from '../../spec/editor.js';
import type {SpecNodeRef} from './node-ref.js';
import {historyKey} from './keys.js';
import {BIND_ONLY, EDITABLE, typeOf} from './prop-model.js';
import {REPOINTING} from '../viewers/viewer-control.js';

/** Read at call time like the kill-walk globals: js-api's `PropertyGrid.update` takes no frame, and
 * the designer's viewer has no table in the shell for the Dart side to find on its own. */
const api = globalThis as {
  grok_Viewer_Get_Look?: (dart: unknown) => unknown;
  grok_PropertyGrid_Update?: (dart: unknown, src: unknown, props: unknown[], table?: unknown) => void;
  grok_Widget_Kill?: (element: Element) => void;
};

const UNAVAILABLE = 'The platform property grid is not available on this client';
/** The platform's color editor pops this into `document.body` and never closes it (core ticket). */
const COLOR_POPUP = '.property-grid-item-editor-color-picker-host';
/** Where a row's focus lands: the editor of a bool/int/choice/string row, or the column combobox. */
const EDITORS = ['input', 'select', '.d4-column-selector'];
/** Registered look props a tag keeps out of the grid: `columnNames` is a write-once alias the
 * platform turns into `filters` at build, so a row over it edits nothing the viewer shows. */
const HIDDEN: Record<string, ReadonlySet<string>> = {'u2-viewer-filters': new Set(['columnNames'])};

type Funnel = (patches: SpecPatch[]) => void;

export interface LookGrid {
  root: HTMLElement;
  /** Re-points the grid when the node was re-built; a no-op while the build is the same. */
  refresh(): void;
}

/** The grid over `x`'s live look, released with `panel`; `funnel` receives what one grid commit
 * becomes — the patches, possibly none. */
export function lookGrid(x: SpecNodeRef, panel: Scope, funnel: Funnel): LookGrid {
  if (api.grok_PropertyGrid_Update === undefined || api.grok_Viewer_Get_Look === undefined)
    return {root: span(UNAVAILABLE, 'u2-designer-hint'), refresh: () => {}};
  const grid = new DG.PropertyGrid();
  grid.root.classList.add('u2-designer-look');
  let current: unknown;
  let attached: Scope | undefined;
  // the textbox editor removes itself on Enter BEFORE the look is written (editor_textbox.dart
  // hideEditor → finishEditing), which drops focus to the body: the committed prop names the row
  // focus goes back to, or the next Ctrl+Z is the platform's ("Nothing to undo")
  let edited: string | undefined;
  const refresh = (): void => {
    const built = x.built();
    if (built === current)
      return;
    attached?.dispose();
    closePopups();
    current = built;
    const wrapper = grid.root.querySelector('.property-grid-wrapper');
    const scrollTop = wrapper?.scrollTop ?? 0;
    const focused = grid.root.contains(document.activeElement) ?
      document.activeElement!.closest<HTMLElement>('tr[data-prop-name]')?.dataset.propName :
      document.activeElement === document.body ? edited : undefined;
    edited = undefined;
    attached = built instanceof DG.Viewer && Component.is(built) ?
      attach(grid, x, built, funnel, (name) => edited = name) : undefined;
    if (attached === undefined)
      grid.update({}, []);
    if (wrapper)
      wrapper.scrollTop = scrollTop;
    for (const row of Array.from(grid.root.querySelectorAll<HTMLElement>('tr[data-prop-name]')))
      row.classList.toggle('u2-designer-look-set', x.node.props?.[row.dataset.propName!] !== undefined);
    if (focused !== undefined)
      focusRow(grid.root.querySelector<HTMLElement>(`tr[data-prop-name="${focused}"]`));
  };
  // the grid's editors are text entry, which the designer's own key map leaves to the browser
  const keydown = (e: Event): void => {
    const run = x.editor === undefined ? undefined : historyKey(e as KeyboardEvent, x.editor);
    if (run === undefined)
      return;
    e.preventDefault();
    e.stopPropagation();
    run();
  };
  grid.root.addEventListener('keydown', keydown);
  panel.own(() => {
    grid.root.removeEventListener('keydown', keydown);
    attached?.dispose();
    closePopups();
    // drops the item controllers' look subscriptions (property_item_view.dart:45-50); the kill
    // releases the Dart widget itself, or every closed panel leaves a PropertyGrid in Widget.getAll()
    grid.update({}, []);
    api.grok_Widget_Kill?.(grid.root);
  });
  refresh();
  return {root: grid.root, refresh};
}

function closePopups(): void {
  for (const popup of Array.from(document.body.querySelectorAll(COLOR_POPUP)))
    popup.remove();
}

/** Focus onto a rebuilt row: its editor where the row shows one, the row itself otherwise (a
 * string row shows a label until clicked) — inside the grid root either way, where the undo chord
 * is the designer's. */
function focusRow(row: HTMLElement | null): void {
  if (row === null)
    return;
  const editor = row.querySelector<HTMLElement>(EDITORS.join(', '));
  if (editor === null)
    row.tabIndex = -1;
  (editor ?? row).focus();
}

/** The registered look props, in registry order — what the document accepts, so what the grid
 * offers; never a bind-only frame (the platform's own `table` row), nor what the tag hides. */
function shownProps(x: SpecNodeRef, viewer: DG.Viewer): DG.Property[] {
  const hidden = HIDDEN[x.node.tag];
  const names = new Set((x.meta()?.props ?? [])
    .filter((p) => typeOf(p) !== BIND_ONLY && hidden?.has(p.name) !== true).map((p) => p.name));
  return viewer.getProperties().filter((p) => names.has(p.name));
}

function update(grid: DG.PropertyGrid, shown: DG.Property[], viewer: DG.Viewer): void {
  const look = api.grok_Viewer_Get_Look!(viewer.dart);
  api.grok_PropertyGrid_Update!(grid.dart, look, shown.map((p) => DG.toDart(p)), viewer.dataFrame?.dart ?? null);
}

/** The grid over `viewer`'s look, and the patch capture: one event per commit (Q3), coalesced per
 * microtask so a multi-property commit is one undo entry; nothing in Run mode (R-c — the document
 * is the authority, and Run-mode edits never reach it). A viewer mid-repoint echoes its
 * construction-time look back, normalized — not an edit; nor is the look's own value at attach
 * time written back (opening the color editor does that). */
function attach(grid: DG.PropertyGrid, x: SpecNodeRef, viewer: DG.Viewer & Component, funnel: Funnel,
  onEdit: (name: string) => void): Scope {
  const scope = new Scope();
  const shown = shownProps(x, viewer);
  update(grid, shown, viewer);
  const look = api.grok_Viewer_Get_Look!(viewer.dart);
  const baseline = new Map(shown.map((p) => [p.name, p.get(look)]));
  const index = new Map(viewer.getProperties().map((p) => [p.name, p]));
  const pending = new Map<string, unknown>();
  const subscription = viewer.onPropertyValueChanged.subscribe((e) => {
    const name = (e as {args?: {property?: {name?: string}}})?.args?.property?.name;
    if (name === undefined || REPOINTING.has(viewer) || !x.instance.designTime)
      return;
    onEdit(name);
    const first = pending.size === 0;
    // the value at event time, off the look as it is now — a later apply re-creates the viewer
    pending.set(name, index.get(name)?.get(api.grok_Viewer_Get_Look!(viewer.dart)));
    if (first)
      queueMicrotask(() => funnel(patchesOf(x, pending, baseline)));
  });
  scope.own(() => subscription.unsubscribe());
  return scope;
}

/** What the grid edits become: nothing for a bound prop (the binding is what the document says),
 * nothing for a prop the document cannot carry (an object), nothing that is not a change — from
 * the document's value or from the look's at attach; `null` — a picker cleared — is "unset", the
 * key removed. */
function patchesOf(x: SpecNodeRef, pending: Map<string, unknown>, baseline: Map<string, unknown>): SpecPatch[] {
  const node = x.node;
  const patches: SpecPatch[] = [];
  const meta = x.meta();
  for (const [name, raw] of pending) {
    const prop = meta?.props.find((p) => p.name === name);
    if (!prop || !EDITABLE.has(typeOf(prop)) || node.bind?.[name] !== undefined)
      continue;
    if (baseline.has(name) && Component.sameValue(raw, baseline.get(name)))
      continue;
    const value = raw === null ? undefined : raw;
    const current = node.props?.[name];
    if (Component.sameValue(current, value) || (current === undefined && value === ''))
      continue;
    patches.push({op: 'set-prop', node, name, value});
  }
  pending.clear();
  return patches;
}
