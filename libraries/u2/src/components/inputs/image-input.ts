import {Input, InputOptions} from '../../core/input-base.js';

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

  get alt(): string {
    return this.options.alt ?? '';
  }

  set alt(x: string) {
    this.options.alt = x ?? undefined;
    this.root.querySelector<HTMLImageElement>('.u2-image-preview')!.alt = this.alt;
  }

  protected createEditor(): HTMLElement {
    const image = document.createElement('img');
    image.className = 'u2-image-preview';
    image.alt = this.options.alt ?? '';

    const remove = document.createElement('button');
    remove.type = 'button';
    remove.className = 'u2-image-remove';
    remove.textContent = 'REMOVE';
    this.scope.listen(remove, 'click', () => this.value = '');

    this.effect(() => {
      const value = this.value;
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
