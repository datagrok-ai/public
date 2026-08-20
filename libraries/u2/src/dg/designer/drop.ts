/* The platform's drag channel, pointed at the tray: a file or a function dragged out of Browse
   becomes a data source. The reading half is pure — the platform objects are read once, into plain
   records — so the payload→patch mapping is testable without a drag. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {nameForTag} from '../../spec/editor.js';
import type {Registry} from '../../spec/registry.js';
import type {SpecNode} from '../../spec/spec.js';
import {sourceNode, specName} from './tray.js';

/** The core function that reads a server file as a table: `fullPath`, `sheetName`, and one
 * `dataframe` output (xamgle `commands.xcmd.g.dart:174`). */
export const OPEN_FILE = 'OpenFile';

const FUNC_SOURCE = 'u2-func-source';
/** What the platform opens as a table with no importer; everything else comes from a file handler. */
const TEXT_TABLE = ['csv', 'tsv', 'txt', 'd42'];

/** One dropped thing, reduced to what a `u2-func-source` needs. */
export interface DropItem {
  kind: 'file' | 'func';
  /** What the node's `func` prop says. */
  ref: string;
  /** What the node's `params` prop says — the path for a file, nothing for a function. */
  params: Record<string, unknown>;
  /** What the component is named after, before uniquification. */
  label: string;
  /** Whether the drop may run it once: a file read is known to be a read (DD9). */
  runs: boolean;
}

export interface DropReading {
  items: DropItem[];
  /** Recognized but unusable, each with its reason — one message per drop. */
  refused: string[];
}

/** Reads the platform's drag payload — one object, or the array a multi-selection drags
 * (`data_source_card_view.dart:153`). Anything that is neither a file nor a function is ignored. */
export function readDrop(dragObject: unknown): DropReading {
  const items: DropItem[] = [];
  const refused: string[] = [];
  let tabular: Set<string> | undefined;
  for (const one of Array.isArray(dragObject) ? dragObject : [dragObject]) {
    if (one instanceof DG.FileInfo) {
      tabular ??= tabularExtensions();
      const refusal = fileRefusal(one, tabular);
      if (refusal === null)
        items.push(fileItem(one));
      else
        refused.push(refusal);
    } else if (one instanceof DG.Func)
      items.push(funcItem(one));
  }
  return {items, refused};
}

/** Extensions the platform can open as a table: the text ones, `d42`, and whatever the registered
 * file handlers declare. Computed per call — a package staged later in the session counts. */
export function tabularExtensions(): Set<string> {
  const exts = new Set(TEXT_TABLE);
  // `file-handler` is matched as tag OR camelCased role, which is the form packages register
  for (const handler of DG.Func.find({meta: {role: 'file-handler'}})) {
    for (const ext of String(handler.options['ext'] ?? '').split(',')) {
      const one = ext.trim().toLowerCase();
      if (one !== '')
        exts.add(one);
    }
  }
  return exts;
}

/** What a spec writes for a function: the bare name, qualified `Namespace:Name` only where the bare
 * one would name more than one function (the WO-16 rule). `nqName` already carries the colon. */
export function funcRef(func: DG.Func): string {
  return DG.Func.find({name: func.name}).length > 1 ? func.nqName : func.name;
}

/** The component node an item becomes. Pure; `unique` reserves the name across both trees. */
export function dropNode(item: DropItem, registry: Registry,
  unique: (base: string) => string): SpecNode {
  const props: Record<string, unknown> = {func: item.ref};
  if (Object.keys(item.params).length > 0)
    props.params = item.params;
  return sourceNode(registry.get(FUNC_SOURCE), FUNC_SOURCE,
    unique(specName(item.label, nameForTag(FUNC_SOURCE))), props);
}

export interface DesignerDropOptions {
  element: HTMLElement;
  /** Whether a drag may land right now: design mode, not disposed, no designer drag in flight. */
  active: () => boolean;
  /** On for the whole drag, so the tray can say where the drop will land. */
  onDragActive: (on: boolean) => void;
  onDrop: (reading: DropReading) => void;
}

/** Wires `element` onto the platform drag channel. `dropIndication: false` is mandatory: the
 * default overlay is body-level, sized to the element's bounding rect, and owns the only mouseup
 * listener (d4 `ui.dart:688-723`). The registration is dart-side and has no counterpart to undo
 * (`ui.dart:737`), so it outlives the designer — every callback asks `active()` first. */
export function makeDesignerDroppable(options: DesignerDropOptions): void {
  ui.makeDroppable(options.element, {
    dropIndication: false,
    acceptDrag: (args: DG.DragDropArgs<unknown>) =>
      options.active() && recognizes(readDrop(args.dragObject)),
    onBeginDrag: () => options.onDragActive(true),
    onEndDrag: () => options.onDragActive(false),
    doDrop: (args: DG.DragDropArgs<unknown>) => {
      options.onDragActive(false);
      if (!options.active())
        return;
      // the platform's own veto flag; never stopPropagation — the drag terminator is a bubbling
      // document mouseup (`html_utils.dart:319`) and the ghost would never come down
      args.handled = true;
      options.onDrop(readDrop(args.dragObject));
    },
  });
}

/** A drag is ours when it carries anything we recognize, a file we will refuse included: an
 * `acceptDrag` that says no makes the drop fall through in silence, which is the worst possible
 * answer to "why did nothing happen". */
function recognizes(reading: DropReading): boolean {
  return reading.items.length > 0 || reading.refused.length > 0;
}

/** Read field by field, never spread: a platform entity keeps everything on prototype getters. */
function fileItem(file: DG.FileInfo): DropItem {
  const name = file.fileName;
  const ext = file.extension;
  const cut = ext === '' ? name.length : name.length - ext.length - 1;
  return {kind: 'file', ref: OPEN_FILE, params: {fullPath: file.fullPath},
    label: name.slice(0, cut), runs: true};
}

function funcItem(func: DG.Func): DropItem {
  return {kind: 'func', ref: funcRef(func), params: {},
    label: func.friendlyName ?? func.name, runs: false};
}

function fileRefusal(file: DG.FileInfo, tabular: Set<string>): string | null {
  const name = file.fileName;
  if (file.isDirectory)
    return `${name}: a folder is not a table — drop a file`;
  const ext = file.extension.toLowerCase();
  if (ext === '')
    return `${name}: nothing opens a file with no extension as a table`;
  return tabular.has(ext) ? null : `${name}: nothing opens .${ext} files as a table`;
}
