import {Input, InputOptions} from '../../core/input-base.js';
import {TypeAhead} from '../inputs/typeahead.js';
import {valueEquals} from '../../core/filter/index.js';
import type {FilterRef} from '../../core/filter/model.js';
import type {FilterProperty, FilterSchema, FilterValueItem} from '../../core/filter/schema.js';

export interface RefInputOptions extends InputOptions<FilterRef | null> {
  prop: FilterProperty;
  /** Its `values` is the type-ahead source, its `renderer` the row presentation. */
  schema: FilterSchema;
  placeholder?: string;
  debounceMs?: number;
}

/** A ref value picked by type-ahead over the schema's value items; the value is the
 * `FilterRef`, the text its display name. */
export class RefInput extends Input<FilterRef | null, RefInputOptions> {
  constructor(options: RefInputOptions) {
    super(options, null);
    this.root.dataset.u2 = 'ref-input';
  }

  protected createEditor(): HTMLElement {
    const {prop, schema} = this.options;
    const renderer = schema.renderer?.(prop);
    const item = (value: FilterRef | null): FilterValueItem | null =>
      value === null ? null : {value, label: value.name ?? value.id};
    const typeAhead = new TypeAhead<FilterValueItem>({
      source: (query, signal) => schema.values!(prop, query, signal),
      itemText: (i) => i.label ?? renderer?.caption(i.value) ?? String(i.value),
      render: renderer?.listItem ? (i) => renderer.listItem!(i.value) : undefined,
      placeholder: this.options.placeholder,
      debounceMs: this.options.debounceMs,
    });
    typeAhead.selected.value = item(this.value.peek());
    this.effect(() => {
      const picked = typeAhead.selected.value;
      const value = picked === null ? null : picked.value as FilterRef;
      if (!valueEquals(value ?? undefined, this.value.peek() ?? undefined))
        this.value.value = value;
    });
    this.effect(() => {
      const value = this.value.value;
      const current = typeAhead.selected.peek()?.value as FilterRef | undefined;
      if (!valueEquals(value ?? undefined, current))
        typeAhead.selected.value = item(value);
    });
    return typeAhead.root;
  }
}
