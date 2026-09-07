import {Input, InputOptions} from '../../core/input-base.js';

export interface BoolInputOptions extends InputOptions<boolean> {
  /** Compact toggle instead of a checkbox. */
  switch?: boolean;
}

export class BoolInput extends Input<boolean, BoolInputOptions> {
  constructor(options: BoolInputOptions = {}) {
    super(options, false);
    this.root.dataset.u2 = 'bool-input';
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
    this.effect(() => box.checked = this.value.value);
    this._listen(box, 'change', () => this.value.value = box.checked);
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
      const on = this.value.value;
      track.classList.toggle('u2-input-switch-on', on);
      track.setAttribute('aria-checked', String(on));
    });
    this._listen(track, 'click', () => this._toggle());
    this._listen(track, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key !== ' ' && key !== 'Enter')
        return;
      e.preventDefault();
      this._toggle();
    });
    return track;
  }

  private _toggle(): void {
    if (this.enabled)
      this.value.value = !this.value.peek();
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
