import {Input, InputOptions} from '../core/input-base.js';

export interface DynamicInputOptions<T> extends InputOptions<T | null> {
  /** Draws the current value; null leaves the editor empty. */
  render: (value: T | null) => HTMLElement | null;
  /** Fixed editor size in pixels (Dart `DynamicInput.size`). */
  size?: {width?: number, height?: number};
}

/** An input that shows its value instead of editing it — the port of Dart's `DynamicInput`
 * (`dynamic_input.dart:3-41`). Every write swaps the rendered element for a fresh one from
 * `render`; there is no user edit to report, so nothing ever fires from the editor itself and the
 * value only ever changes from code (Dart's `onInput` is an empty stream, `:26`).
 *
 * `render` output is foreign DOM: it is removed on the next write and on dispose, never
 * introspected. A renderer that builds a component owns disposing it (`u2/dg`'s `metaInput` does). */
export class DynamicInput<T = unknown> extends Input<T | null, DynamicInputOptions<T>> {
  private _box!: HTMLElement;
  // assigned in createEditor, which the base constructor runs before field initializers
  private _rendered!: HTMLElement | null;

  constructor(options: DynamicInputOptions<T>) {
    super(options, null);
    this.root.dataset.u2 = 'dynamic-input';
  }

  protected createEditor(): HTMLElement {
    this._rendered = null;
    const editor = document.createElement('div');
    this._box = editor;
    editor.className = 'u2-dynamic';
    const size = this.options.size;
    if (size?.width !== undefined)
      editor.style.width = `${size.width}px`;
    if (size?.height !== undefined)
      editor.style.height = `${size.height}px`;

    this.effect(() => this._show(this.options.render(this.value.value)));
    this.own(() => this._show(null));
    return editor;
  }

  private _show(el: HTMLElement | null): void {
    this._rendered?.remove();
    this._rendered = el;
    if (el !== null)
      this._box.append(el);
  }
}
