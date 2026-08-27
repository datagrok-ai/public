/** KNIME-style colored regions behind the graph — purely visual, owned by FlowEditor,
 *  rendered inside `area.content.holder` so pan/zoom move them with the graph. */

const ANN_DEFAULT_COLOR = '#BBDEFB';
const ANN_DEFAULT_BORDER = '#1976D2';
export const ANN_DEFAULT_FONT_SIZE = 13;

export interface AnnotationDoc {
  id: string;
  pos: {x: number; y: number};
  size: {w: number; h: number};
  text: string;
  color: string;
  /** Title font size in px; absent = the default (13). */
  fontSize?: number;
}

/** Material 100-level backgrounds paired with their 700-level borders. */
export const ANNOTATION_COLORS: Array<{name: string; bg: string; border: string}> = [
  {name: 'Blue', bg: '#BBDEFB', border: '#1976D2'},
  {name: 'Indigo', bg: '#C5CAE9', border: '#303F9F'},
  {name: 'Cyan', bg: '#B2EBF2', border: '#0097A7'},
  {name: 'Teal', bg: '#B2DFDB', border: '#00796B'},
  {name: 'Green', bg: '#C8E6C9', border: '#388E3C'},
  {name: 'Lime', bg: '#F0F4C3', border: '#9E9D24'},
  {name: 'Yellow', bg: '#FFF59D', border: '#F9A825'},
  {name: 'Orange', bg: '#FFE0B2', border: '#EF6C00'},
  {name: 'Red', bg: '#FFCDD2', border: '#D32F2F'},
  {name: 'Pink', bg: '#F8BBD0', border: '#C2185B'},
  {name: 'Purple', bg: '#E1BEE7', border: '#7B1FA2'},
  {name: 'Brown', bg: '#D7CCC8', border: '#5D4037'},
  {name: 'Gray', bg: '#ECEFF1', border: '#607D8B'},
];

export const ANNOTATION_TITLE_SIZES: Array<{name: string; size: number}> = [
  {name: 'Small', size: 11},
  {name: 'Normal', size: ANN_DEFAULT_FONT_SIZE},
  {name: 'Medium', size: 16},
  {name: 'Large', size: 20},
  {name: 'Huge', size: 26},
];

function borderForBackground(bg: string): string {
  const match = ANNOTATION_COLORS.find((c) => c.bg.toLowerCase() === bg.toLowerCase());
  return match ? match.border : ANN_DEFAULT_BORDER;
}

export class FlowAnnotation {
  id: string;
  pos: {x: number; y: number};
  size: {w: number; h: number};
  text: string;
  color: string;
  fontSize: number;
  /** Outer wrapper element — added to `area.content.holder` by the editor. */
  readonly element: HTMLElement;
  readonly titleEl: HTMLElement;
  readonly resizeHandle: HTMLElement;

  constructor(opts: Partial<AnnotationDoc> = {}) {
    this.id = opts.id ?? `ann-${Math.random().toString(36).slice(2, 10)}`;
    this.pos = opts.pos ?? {x: 0, y: 0};
    this.size = opts.size ?? {w: 240, h: 140};
    this.text = opts.text ?? 'Annotation';
    this.color = opts.color ?? ANN_DEFAULT_COLOR;
    const fs = Number(opts.fontSize);
    this.fontSize = Number.isFinite(fs) && fs > 0 ? fs : ANN_DEFAULT_FONT_SIZE;

    this.element = document.createElement('div');
    this.element.className = 'ff-annotation';
    this.element.dataset.annotationId = this.id;

    this.titleEl = document.createElement('div');
    this.titleEl.className = 'ff-annotation-title';
    this.titleEl.textContent = this.text;
    this.titleEl.contentEditable = 'true';
    this.titleEl.spellcheck = false;
    this.titleEl.addEventListener('input', () => {
      this.text = this.titleEl.textContent ?? '';
    });
    this.element.appendChild(this.titleEl);

    this.resizeHandle = document.createElement('div');
    this.resizeHandle.className = 'ff-annotation-resize';
    this.element.appendChild(this.resizeHandle);

    this.applyAll();
  }

  applyAll(): void {
    this.applyPos();
    this.applySize();
    this.applyColor();
    this.applyFont();
  }
  applyPos(): void {
    this.element.style.left = `${this.pos.x}px`;
    this.element.style.top = `${this.pos.y}px`;
  }
  applySize(): void {
    this.element.style.width = `${this.size.w}px`;
    this.element.style.height = `${this.size.h}px`;
  }
  applyColor(): void {
    this.element.style.background = this.color;
    this.element.style.borderColor = borderForBackground(this.color);
  }
  applyFont(): void {
    this.titleEl.style.fontSize = `${this.fontSize}px`;
  }

  toDoc(): AnnotationDoc {
    return {
      id: this.id,
      pos: {...this.pos},
      size: {...this.size},
      text: this.text,
      color: this.color,
      // Omitted at the default so untouched saves stay tidy (same pattern as autorun).
      ...(this.fontSize !== ANN_DEFAULT_FONT_SIZE ? {fontSize: this.fontSize} : {}),
    };
  }
}
