/* `domains.search` — the ribbon search box over a source: the search variant of `TextInput`
   (magnifier, ✕), writing `source.search` after a quiet spell, on Enter at once, cleared by
   Escape; two-way, so a search set in code shows in the box. The server runs it over the
   table's `searchable` columns, AND-ed with the query. A search the USER makes goes through the
   gate (STATE-CONTRACT H6), as the filters do — it answers at once while there is nothing to lose
   and refuses while a batch is being written back; cancel puts the previous text back. */
import {computed, TextInput} from '@datagrok-libraries/u2';
import type {DomainSource} from '../core/domain-source.js';
import {confirmDiscard} from '../core/session.js';
import {DomainTable} from './facade.js';

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
        const columns = DomainSearch.columns(source).join(', ');
        // a list too long for the box is left to the tooltip
        return columns === '' || columns.length > 20 ? `Search ${DomainTable.noun(source.schema.info)}…` :
          `Search ${columns}…`;
      })});
    this.root.classList.add('u2-domain-search');
    this.root.dataset.u2 = 'domain-search';
    const box = this.root.querySelector('input')!;
    this.effect(() => {
      source.state.value;
      const columns = DomainSearch.columns(source);
      box.title = columns.length === 0 ? '' : `Finds the ${DomainTable.noun(source.schema.info)} whose ` +
        `${columns.join(', ')} contain the text — within what the filter shows`;
    });
    const wait = options.debounceMs ?? 300;
    const write = () => {
      this._cancel();
      const text = this.value.peek();
      if (source.search.peek() === text)
        return;
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

  /** The searched columns as the form captions them, in lower case. */
  static columns(source: DomainSource): string[] {
    const props = source.schema.properties;
    return source.schema.info.searchableColumns.map((name) => {
      const prop = props.find((p) => p.name === name);
      return (prop?.friendlyName ?? prop?.caption ?? name.replace(/_/g, ' ')).toLowerCase();
    });
  }

  private _cancel(): void {
    if (this._timer !== undefined)
      clearTimeout(this._timer);
    this._timer = undefined;
  }
}
