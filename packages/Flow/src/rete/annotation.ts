/** Workflow annotations — KNIME-style colored regions sitting *behind* the
 *  graph. Purely visual — they're not part of the Rete editor (no nodes, no
 *  ports), they don't affect compilation, and they live in their own
 *  `Map<id, FlowAnnotation>` on the FlowEditor.
 *
 *  Each annotation owns a `<div>` rendered inside `area.content.holder`, the
 *  transformed canvas content layer, so panning and zooming move them with
 *  the rest of the graph. Stacking is handled in CSS — `z-index: 0` keeps
 *  them under the nodes (which the order extension drives up). */

const ANN_DEFAULT_COLOR = '#BBDEFB';   // pale blue body
const ANN_DEFAULT_BORDER = '#1976D2';  // blue border

export interface AnnotationDoc {
  id: string;
  pos: {x: number; y: number};
  size: {w: number; h: number};
  text: string;
  color: string;
}

export const ANNOTATION_COLORS: Array<{name: string; bg: string; border: string}> = [
  {name: 'Blue',   bg: '#BBDEFB', border: '#1976D2'},
  {name: 'Yellow', bg: '#FFF59D', border: '#F9A825'},
  {name: 'Green',  bg: '#C8E6C9', border: '#388E3C'},
  {name: 'Pink',   bg: '#F8BBD0', border: '#C2185B'},
  {name: 'Gray',   bg: '#ECEFF1', border: '#607D8B'},
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

  toDoc(): AnnotationDoc {
    return {
      id: this.id,
      pos: {...this.pos},
      size: {...this.size},
      text: this.text,
      color: this.color,
    };
  }
}
