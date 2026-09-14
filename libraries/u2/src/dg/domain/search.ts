/* `domains.search` — the ribbon search box over a source: the search variant of `TextInput`
   (magnifier, ✕), writing `source.search` after a quiet spell, on Enter at once, cleared by
   Escape; two-way, so a search set in code shows in the box. The server runs it over the
   table's `searchable` columns, AND-ed with the query. A search the USER makes while the session
   holds unsaved changes goes through the gate (STATE-CONTRACT H6), as the filters do: cancel puts
   the previous text back. */
import {computed} from '../../core/signals.js';
import {TextInput} from '../../components/inputs/text-input.js';
import type {DomainSource} from '../../sources/domain-source.js';
import {confirmDiscard} from '../../sources/session.js';

export interface DomainSearchOptions {
  placeholder?: string;
  /** How long the box waits after the last keystroke before it searches (default 300 ms). */
  debounceMs?: number;
}

export class DomainSearch extends TextInput {
  private _timer: ReturnType<typeof setTimeout> | undefined;

  constructor(readonly source: DomainSource, options: DomainSearchOptions = {}) {
    super({search: true, inline: true, name: 'search', value: source.search.peek(),
      placeholder: options.placeholder ?? computed(() => {
        source.state.value;
        return `Search ${source.schema.info.pluralName.toLowerCase() || 'rows'}…`;
      })});
    this.root.dataset.u2 = 'domain-search';
    const wait = options.debounceMs ?? 300;
    const write = () => {
      this._cancel();
      const text = this.value.peek();
      if (source.search.peek() === text)
        return;
      if (!source.session.isDirty.peek()) {
        source.search.value = text;
        return;
      }
      void confirmDiscard(source.session, {action: 'change the search'}).then((ok) => {
        if (ok)
          source.search.value = text;
        else
          this.value.value = source.search.peek();
      });
    };
    this.effect(() => {
      const text = this.value.value;
      this._cancel();
      if (text !== source.search.peek())
        this._timer = setTimeout(write, wait);
    });
    this.effect(() => {
      const text = source.search.value;
      if (this.value.peek() !== text) {
        this._cancel();
        this.value.value = text;
      }
    });
    const input = this.root.querySelector('input')!;
    const onKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Enter')
        write();
      else if (e.key === 'Escape') {
        this.value.value = '';
        write();
      } else
        return;
      e.preventDefault();
    };
    input.addEventListener('keydown', onKeyDown);
    this.own(() => {
      input.removeEventListener('keydown', onKeyDown);
      this._cancel();
    });
  }

  private _cancel(): void {
    if (this._timer !== undefined)
      clearTimeout(this._timer);
    this._timer = undefined;
  }
}
