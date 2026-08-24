/* The designer's keyboard map: undo/redo, Delete, and Escape — which a drag in flight owns. */
import type {ReadonlySignal} from '../../core/signals.js';
import type {SpecEditor} from '../../spec/editor.js';
import type {SpecInstance, SpecNode} from '../../spec/spec.js';
import {DELETE} from './actions.js';

/** Elements that own the keystroke themselves: the tree's rename box, the Open dialog, the panel. */
const TEXT_ENTRY = /^(INPUT|TEXTAREA|SELECT)$/;

export interface KeyHost {
  readonly editor: ReadonlySignal<SpecEditor | undefined>;
  instance(): SpecInstance | undefined;
  selected(): SpecNode | null;
  multi(): SpecNode[];
  inDrag(): boolean;
  run(name: string): void;
  collapse(): void;
  select(node: SpecNode): void;
  endDrag(commit: boolean): void;
}

export function keyDown(e: KeyboardEvent, host: KeyHost): void {
  const target = e.target as HTMLElement | null;
  if (target && (TEXT_ENTRY.test(target.tagName) || target.isContentEditable))
    return;
  const editor = host.editor.peek();
  if (!editor)
    return;
  const ctrl = e.ctrlKey || e.metaKey;
  let run = historyKey(e, editor);
  if (run === undefined && e.key === 'Delete' && !ctrl)
    run = () => host.run(DELETE);
  else if (run === undefined && e.key === 'Escape' && !host.inDrag()) {
    // a drag in flight owns Escape (the document-level `dragKey` cancels it); a
    // multi-selection collapses to its lead first, the walk to the parent comes after (F5)
    if (host.multi().length > 1)
      run = () => host.collapse();
    else {
      const selected = host.selected();
      const parent = selected ? host.instance()?.parentOf(selected) : null;
      if (parent)
        run = () => host.select(parent);
    }
  }
  if (!run)
    return;
  e.preventDefault();
  e.stopPropagation();
  run();
}

/** The undo/redo chord — Ctrl+Z, Ctrl+Y, Ctrl+Shift+Z — as the call it makes; undefined for any other key. */
export function historyKey(e: KeyboardEvent, editor: SpecEditor): (() => void) | undefined {
  const ctrl = e.ctrlKey || e.metaKey;
  const key = e.key.toLowerCase();
  if (ctrl && key === 'z' && !e.shiftKey)
    return () => editor.undo();
  return ctrl && (key === 'y' || (key === 'z' && e.shiftKey)) ? () => editor.redo() : undefined;
}

export function dragKey(e: KeyboardEvent, host: Pick<KeyHost, 'inDrag' | 'endDrag'>): void {
  if (host.inDrag() && e.key === 'Escape')
    host.endDrag(false);
}
