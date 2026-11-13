import {Input, InputOptions} from '../../core/input-base.js';

export interface BoolInputOptions extends InputOptions<boolean> {
  /** Compact toggle instead of a checkbox. */
  switch?: boolean;
  /** Checklist row: the box first, its caption after it, no label column reserved; a click on the
   * caption toggles. */
  captionAfter?: boolean;
}

export class BoolInput extends Input<boolean, BoolInputOptions> {
  constructor(options: BoolInputOptions = {}) {
    super(options, false);
    this.root.dataset.u2 = 'bool-input';
    if (options.captionAfter)
      this.root.classList.add('u2-input-caption-after');
    // the checkbox's caption is its `<label for>`; the switch is a span a label cannot target
    if (options.switch) {
      this.scope.listen(this.root, 'click', (e) => {
        if ((e.target as Element).closest('[data-u2-part="label"]') !== null)
          this._toggle();
      });
    }
  }

  /** The base state plus `checked`, the value as the box shows it. */
  getState(): ReturnType<Input<boolean>['getState']> {
    return {...super.getState(), checked: this.peek() === true};
  }

  protected createEditor(): HTMLElement {
    const editor = document.createElement('span');
    editor.className = 'u2-input-bool';
    editor.append(this.options.switch ? this._createSwitch() : this._createCheckbox());
    return editor;
  }

  private _createCheckbox(): HTMLElement {
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.className = 'u2-input-checkbox';
    this.effect(() => box.checked = this.value === true);
    this.scope.listen(box, 'change', () => this.value = box.checked);
    return box;
  }

  private _createSwitch(): HTMLElement {
    const track = document.createElement('span');
    track.className = 'u2-input-switch';
    track.tabIndex = 0;
    track.setAttribute('role', 'switch');
    const thumb = document.createElement('span');
    thumb.className = 'u2-input-switch-thumb';
    track.append(thumb);
    this.effect(() => {
      const on = this.value === true;
      track.classList.toggle('u2-input-switch-on', on);
      track.setAttribute('aria-checked', String(on));
    });
    this.scope.listen(track, 'click', () => this._toggle());
    this.scope.listen(track, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key !== ' ' && key !== 'Enter')
        return;
      e.preventDefault();
      this._toggle();
    });
    return track;
  }

  private _toggle(): void {
    if ((this.disabled === null) && !this.readOnly)
      this.value = !this.peek();
  }
}
