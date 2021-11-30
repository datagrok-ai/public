import {ChemPalette} from './chem-palette';
import * as DG from 'datagrok-api/dg';

const cp = new ChemPalette('grok');

export function expandColumn(col:DG.Column,
  grid:DG.Grid,
  cellRenderSize: (celVal:string)=>number,
  textSizeMult = 10,
  minSize = 30,
  maxSize = 650,
  timeout = 500) {
  let maxLen = 0;
  col.categories.forEach((ent: string) => {
    const len = cellRenderSize(ent);
    if (len > maxLen) {
      maxLen = len;
    }
  });
  setTimeout(() => {
      grid.columns.byName(col.name)!.width = Math.min(Math.max(maxLen * textSizeMult, minSize), maxSize);
  },
  timeout);
}

export function setAARRenderer(col: DG.Column, grid: DG.Grid | null = null, grouping = false) {
  col.semType = 'aminoAcids';
  col.setTag('cell.renderer', 'aminoAcids');
  if (grouping) {
    col.setTag('groups', `${grouping}`);
  }
  if (grid) {
    expandColumn(col, grid, (ent) => measureAAR(ent));
  }
}

export function measureAAR(s: string): number {
  const end = s.lastIndexOf(')');
  const beg = s.indexOf('(');
  return end == beg ? s.length : s.length - (end - beg) + 1;
}

function printLeftCentered(
  x: number, y: number, w: number, h: number,
  g: CanvasRenderingContext2D, s: string, color = ChemPalette.undefinedColor,
  pivot: number = 0, left = false, hideMod = false) {
  g.textAlign = 'start';
  let colorPart = pivot == -1 ? s.substring(0) : s.substring(0, pivot);
  if (colorPart.length == 1) {
    colorPart = colorPart.toUpperCase();
  }
  if (colorPart.length >= 3) {
    if (colorPart.substring(0, 3) in ChemPalette.AAFullNames) {
      colorPart = ChemPalette.AAFullNames[s.substring(0, 3)] + colorPart.substr(3);
    } else if (colorPart.substring(1, 4) in ChemPalette.AAFullNames) {
      colorPart = colorPart[0] + ChemPalette.AAFullNames[s.substring(1, 4)] + colorPart.substr(4);
    }
  }
  let grayPart = pivot == -1 ? '' : s.substr(pivot);
  if (hideMod) {
    let end = colorPart.lastIndexOf(')');
    let beg = colorPart.indexOf('(');
    if (beg > -1 && end > -1 && end - beg > 2) {
      colorPart = colorPart.substr(0, beg) + '(+)' + colorPart.substr(end + 1);
    }

    end = grayPart.lastIndexOf(')');
    beg = grayPart.indexOf('(');
    if (beg > -1 && end > -1 && end - beg > 2) {
      grayPart = grayPart.substr(0, beg) + '(+)' + grayPart.substr(end + 1);
    }
  }
  const textSize = g.measureText(colorPart + grayPart);
  const indent = 5;

  const colorTextSize = g.measureText(colorPart);

  function draw(dx1: number, dx2: number) {
    g.fillStyle = color;
    g.fillText(
      colorPart,
      x + dx1,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    g.fillStyle = ChemPalette.undefinedColor;
    g.fillText(
      grayPart,
      x + dx2,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
  }


  if (left || textSize.width > w) {
    draw(indent, indent + colorTextSize.width);
    return x + colorTextSize.width + g.measureText(grayPart).width;
  } else {
    draw((w - textSize.width) / 2, (w - textSize.width) / 2 + colorTextSize.width);
    return x + (w - textSize.width) / 2 + colorTextSize.width;
  }
}


export class AminoAcidsCellRenderer extends DG.GridCellRenderer {
    chemPalette: ChemPalette | null;

    get name() {
      return 'aminoAcidsCR';
    }

    get cellType() {
      return 'aminoAcids';
    }

    constructor() {
      super();
      this.chemPalette = null;
    }

    render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
      gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
      if (this.chemPalette === null) {
        this.chemPalette = new ChemPalette('grok', gridCell.tableColumn?.getTag('groups') ? true : false);
      }
      g.save();
      g.beginPath();
      g.rect(x, y, w, h);
      g.clip();
      g.font = `14px monospace`;
      g.textBaseline = 'top';
      const s: string = gridCell.cell.value ? gridCell.cell.value : '-';
      const [color, pivot] = cp.getColorPivot(s);
      printLeftCentered(x, y, w, h, g, s, color, pivot, false, true);
      g.restore();
    }
}


export class AlignedSequenceCellRenderer extends DG.GridCellRenderer {
  get name() {
    return 'alignedSequenceCR';
  }

  get cellType() {
    return 'alignedSequence';
  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle ) {
    w = Math.min(gridCell.grid.canvas.width - x, w);
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    g.font = '14px monospace';
    g.textBaseline = 'top';
    const s: string = gridCell.cell.value ?? '';

    //TODO: can this be replaced/merged with splitSequence?
    const subParts = s.split('-');
    const [text, simplified] = processSequence(subParts);
    const textSize = g.measureText(text.join(''));
    x = Math.max(x, x + (w - textSize.width) / 2);

    subParts.forEach((amino: string, index) => {
      const [color, pivot] = cp.getColorPivot(amino);
      g.fillStyle = ChemPalette.undefinedColor;
      if (index + 1 < subParts.length) {
        const gap = simplified?'':' ';
        amino += `${amino?'':'-'}${gap}`;
      }
      x = printLeftCentered(x, y, w, h, g, amino, color, pivot, true);
    });

    g.restore();
  }
}


export function processSequence(subParts:string[]) : [string[], boolean] {
  const simplified = !subParts.some((amino, index) =>
    amino.length > 1 &&
    index != 0 &&
    index != subParts.length - 1);

  const text:string[] = [];
  subParts.forEach((amino: string, index) => {
    if (index < subParts.length) {
      const gap = simplified?'':' ';
      amino += `${amino?'':'-'}${gap}`;
    }
    text.push(amino);
  });
  return [text, simplified];
}
