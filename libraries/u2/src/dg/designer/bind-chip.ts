/* The binding chip (UB-7b): what a bound property row renders in place of its editor — the path,
   `⇄` where the declaration says two-way, `×` to unbind, a click to re-pick. An Input<string>, so
   it rides the form's label column, `data-u2-prop` stamping, `input(name)` lookup and `refresh()`
   for free; `value` holds the path and is the form's to refresh. Commits stay with the caller. */
import {Input, InputOptions} from '../../core/input-base.js';
import {span} from '../../core/elements.js';
import {iconButton} from '../../components/actions/buttons.js';

export interface BindChipOptions extends InputOptions<string> {
  path: string;
  /** Declaration-based (`prop.twoWay === true`) — the chip describes the document, never probes. */
  twoWay: boolean;
  onPick: () => void;
  onClear: () => void;
}

export class BindChip extends Input<string, BindChipOptions> {
  private _clear!: HTMLElement;

  constructor(options: BindChipOptions) {
    super({...options, value: options.path}, '');
    this.root.dataset.u2 = 'bind-chip';
  }

  /** The × button — the caller stamps `dataset.u2Unbind` on it. */
  get clearButton(): HTMLElement {
    return this._clear;
  }

  protected createEditor(): HTMLElement {
    const chip = document.createElement('div');
    chip.className = 'u2-bind-chip';
    const path = span('', 'u2-bind-chip-path');
    chip.append(path);
    if (this.options.twoWay)
      chip.append(span('⇄', 'u2-bind-chip-two-way'));
    this._clear = iconButton('times', this.options.onClear, {tooltip: 'Remove binding'});
    this._clear.classList.add('u2-bind-chip-clear');
    chip.append(this._clear);
    const onClick = (e: Event) => {
      if (!this._clear.contains(e.target as Node))
        this.options.onPick();
    };
    chip.addEventListener('click', onClick);
    this.own(() => chip.removeEventListener('click', onClick));
    // the row's tooltip carries the full path — the chip may ellipsize it
    this.effect(() => this.root.title = path.textContent = this.value.value ?? '');
    return chip;
  }
}
