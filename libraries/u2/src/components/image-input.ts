import {Input, InputOptions} from '../core/input-base.js';

export interface ImageInputOptions extends InputOptions<string> {
  alt?: string;
}

/** Image by URL or data URL, the port of Dart's `ImageInput`: the value is the source string, the
 * editor is the preview, and a REMOVE link revealed on hover clears it. */
export class ImageInput extends Input<string, ImageInputOptions> {
  constructor(options: ImageInputOptions = {}) {
    super(options, '');
    this.root.dataset.u2 = 'image-input';
  }

  protected createEditor(): HTMLElement {
    const image = document.createElement('img');
    image.className = 'u2-image-preview';
    image.alt = this.options.alt ?? '';

    const remove = document.createElement('button');
    remove.type = 'button';
    remove.className = 'u2-image-remove';
    remove.textContent = 'REMOVE';
    const onRemove = () => this.value.value = '';
    remove.addEventListener('click', onRemove);
    this.own(() => remove.removeEventListener('click', onRemove));

    this.effect(() => {
      const value = this.value.value;
      image.src = value;
      image.hidden = value === '';
      remove.hidden = value === '';
    });

    const editor = document.createElement('div');
    editor.className = 'u2-image';
    editor.append(image, remove);
    return editor;
  }
}
