/* The one binding between a collection control and its source (Astra C3): the control's lead row is
   the source's `currentRow`, and the control's multi-selection is the frame's selection bitset —
   gesture-armed, so the lead a load puts on the first row is not a user's selection. `DomainList`
   and `DomainDataTable` share it, and a kanban or a timeline consumes the same contract instead of
   growing a third copy. */
import {signal, untracked} from '../../core/signals.js';
import type {Control} from '../../core/component.js';
import type {VirtualRows} from '../../components/collections/list.js';
import type {DomainSource} from '../../sources/domain-source.js';
import type {RowView} from '../../sources/rows-like.js';

export class DomainSelection {
  /** Binds `rows` (any `VirtualRows`) to `source` inside `owner`'s scope, and answers nothing: the
   * binding lives as long as the owner does.
   *
   * Three rules, each of them what the two copies already did:
   * 1. the lead and `source.currentRow` are one value, written in either direction, seeded from the
   *    row that is current already so the first effect cannot clear it;
   * 2. the multi-selection is mirrored into `df.selection` only after a click or a keypress
   *    (`armed`), and a collection read again keeps it BY KEY — a query that answers none of the
   *    kept keys is a new collection with nothing selected;
   * 3. Escape clears both the mark and the lead. */
  static bind(owner: Control, source: DomainSource, rows: VirtualRows<RowView>): void {
    // the list's selection and the source's current row are one thing, written in either direction —
    // seeded from the row that is current already, so the first effect does not clear it
    const current = source.currentRow.peek();
    if (current !== null)
      rows.selectedIndex.value = source.rows.items.peek().findIndex((r) => r.id === current.id);
    owner.effect(() => {
      const row = source.rows.items.peek()[rows.selectedIndex.value] ?? null;
      if (row?.id !== source.currentRow.peek()?.id)
        source.currentRow.value = row;
    });
    owner.effect(() => {
      const row = source.currentRow.value;
      const items = source.rows.items.value;
      const at = row === null ? -1 : items.findIndex((r) => r.id === row.id);
      if (at !== rows.selectedIndex.peek())
        rows.selectedIndex.value = at;
    });
    // the control's multi-selection IS the collection's: everything that acts on "the selected rows"
    // — `domains.bulkEdit`, `source.stageRestoreSelection` — reads `source.selection`, which follows
    // the frame's selection bitset, and nothing else on a list page writes it. A selection is a
    // USER's: the lead a load puts on the first row is not one, so nothing is mirrored until the
    // control is clicked or keyed, and a new collection starts over.
    const picked = signal(false);
    const onPick = () => picked.value = true;
    const onEscape = (e: Event) => {
      if ((e as KeyboardEvent).key !== 'Escape')
        return;
      picked.value = false;
      rows.selectedIndex.value = -1;
    };
    rows.root.addEventListener('click', onPick);
    rows.root.addEventListener('keydown', onPick);
    rows.root.addEventListener('keydown', onEscape);
    owner.own(() => {
      rows.root.removeEventListener('click', onPick);
      rows.root.removeEventListener('keydown', onPick);
      rows.root.removeEventListener('keydown', onEscape);
    });
    let frame: unknown;
    let kept: string[] = [];
    owner.effect(() => {
      const selected = rows.selectedIndices.value;
      const df = source.df.value;
      const own = picked.value;
      const bits = df?.selection as
        {get?(i: number): boolean, set?(i: number, value: boolean): void} | null | undefined;
      if (df === undefined || typeof bits?.set !== 'function')
        return;
      // a collection read again is the same collection: a selection survives it by KEY (a bulk
      // edit refreshes, and the rows it wrote are still the rows the user picked), and a query
      // that answers none of them is a new collection with nothing selected
      if (df !== frame) {
        frame = df;
        // nothing to carry over leaves the control's own lead alone: a load puts it on the first row
        const again = own ? kept : [];
        if (again.length === 0) {
          picked.value = false;
          return;
        }
        // the scroller copies the items in an effect of its own, which runs AFTER this one: the
        // restore waits for it, or it would look the keys up in the collection that just left
        queueMicrotask(() => {
          if (owner.scope.isDisposed)
            return;
          const present = new Set(source.rows.items.peek().map((row) => row.id));
          const back = again.filter((id) => present.has(id));
          picked.value = back.length > 0;
          if (back.length > 0)
            untracked(() => rows.selectKeys(back));
        });
        return;
      }
      // the emptying the scroller does when its items are replaced is not the user clearing
      // the selection: only a non-empty one is remembered, and `own` is what clears it
      if (own && selected.size > 0)
        kept = rows.selectedKeys();
      // read before writing: a bitset that fires per write would send a change event per ROW on
      // every pass, and a frame the platform writer holds is not to be touched for nothing
      for (let i = 0; i < df.rowCount; i++) {
        const on = own && selected.has(i);
        if (typeof bits.get !== 'function' || bits.get(i) !== on)
          bits.set(i, on);
      }
    });
  }
}
