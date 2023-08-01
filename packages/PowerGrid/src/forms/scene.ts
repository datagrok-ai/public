import * as DG from "datagrok-api/dg";

export interface IStyle {
  color?: string;
  backColor?: string;
  horzAlign?: 'left' | 'right' | 'center';
  vertAlign?: 'top' | 'bottom' | 'center';
}


export abstract class Element {
  bounds: DG.Rect;
  parent?: Element;
  style?: IStyle;

  protected constructor(bounds: DG.Rect) {
    this.bounds = bounds;
  }

  abstract render(g: CanvasRenderingContext2D): void;
}


export class LabelElement extends Element {
  text: string;

  constructor(bounds: DG.Rect, text: string, style?: IStyle) {
    super(bounds);
    this.text = text;
    this.style = style;
  }

  render(g: CanvasRenderingContext2D) {
    g.fillStyle = this.style?.color ?? 'grey'
    g.textAlign = this.style?.horzAlign ?? 'left';
    g.textBaseline = 'middle';
    //g.fillText(this.text, this.bounds.x, this.bounds.y, this.bounds.width);
    drawClipped(g, this.bounds, () => {
      const x
        = g.textAlign == 'left' ? this.bounds.x
        : g.textAlign == 'right' ? this.bounds.right
        : g.textAlign == 'center' ? this.bounds.midX : 0;
      g.fillText(this.text, this.bounds.x, this.bounds.midY, this.bounds.width);
    });
  }
}


export class GridCellElement extends Element {
  gridCell: DG.GridCell;

  constructor(bounds: DG.Rect, gridCell: DG.GridCell) {
    super(bounds);
    this.gridCell = gridCell;
  }

  render(g: CanvasRenderingContext2D) {
    g.fillStyle = 'grey';
    this.gridCell.renderer.render(g,
      this.bounds.x, this.bounds.y, this.bounds.width, this.bounds.height,
      this.gridCell, this.gridCell.style);
  }
}


export class Scene extends Element {
  elements: Element[] = [];

  constructor(bounds: DG.Rect) {
    super(bounds);
  }

  render(g: CanvasRenderingContext2D) {
    for (const e of this.elements) {
      e.render(g);
    }
  }
}


function drawClipped(g: CanvasRenderingContext2D, bounds: DG.Rect, draw: () => void) {
  try {
    g.save();
    g.rect(bounds.x, bounds.y, bounds.width, bounds.height);
    g.clip();
    draw();
  }
  finally {
    g.restore();
  }
}