import {ChemPalette} from './chem-palette';
import * as DG from 'datagrok-api/dg';

const cp =new ChemPalette('grok');


function printLeftCentered(
  x: number,
  y: number,
  w: number,
  h: number,
  g: CanvasRenderingContext2D,
  s: string,
  color = ChemPalette.undefinedColor,
  pivot: number = 0,
  left = false,
  hideMod = false,
) {
  let colorPart = pivot == -1 ? s.substring(0) : s.substring(0, pivot);
  if (colorPart.length == 1) {
    colorPart = colorPart.toUpperCase();
  }
  if (colorPart.length >= 3) {
    if (colorPart.substring(0, 3) in ChemPalette.AAFullNames) {
      colorPart = ChemPalette.AAFullNames[s.substring(0, 3)] + colorPart.substr(3);
    } else if (colorPart.substring(1, 4) in ChemPalette.AAFullNames) {
      colorPart = colorPart.at(0) + ChemPalette.AAFullNames[s.substring(1, 4)] + colorPart.substr(4);
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
  if (left || textSize.width > w) {
    g.fillStyle = color;
    g.fillText(
      colorPart,
      x + indent,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    g.fillStyle = ChemPalette.undefinedColor;
    g.fillText(
      grayPart,
      x + indent + colorTextSize.width,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    return x + colorTextSize.width + g.measureText(grayPart).width;
    ;
  } else {
    g.fillStyle = color;
    g.fillText(
      colorPart,
      x + (w - textSize.width) / 2,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    g.fillStyle = ChemPalette.undefinedColor;
    g.fillText(
      grayPart,
      x + (w - textSize.width) / 2 + colorTextSize.width,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    return x + (w - textSize.width) / 2 + colorTextSize.width;
  }
}

export class AminoAcidsCellRenderer extends DG.GridCellRenderer {
    private fontSize = 15;


    get name() {
      return 'aminoAcidsCR';
    }

    get cellType() {
      return 'aminoAcids';
    }

    get defaultWidth() {
      return 30;
    }


    render(
      g: CanvasRenderingContext2D,
      x: number,
      y: number,
      w: number,
      h: number,
      gridCell: DG.GridCell,
      cellStyle: DG.GridCellStyle,
    ) {
      g.save();
      g.beginPath();
      g.rect(x, y, w, h);
      g.clip();
      g.font = `${this.fontSize}px monospace`;
      g.textBaseline = 'top';
      const s: string = gridCell.cell.value ? gridCell.cell.value : '-';
      const [color, pivot] = cp.getColorPivot(s);
      printLeftCentered(x, y, w, h, g, s, color, pivot, false, true);
      g.restore();
    }
}

export class AlignedSequenceCellRenderer extends DG.GridCellRenderer {
    private maxCellWidth = 270;


    constructor() {
      super();
    }


    get name() {
      return 'alignedSequenceCR';
    }

    get cellType() {
      return 'alignedSequence';
    }

    get defaultWidth() {
      return this.maxCellWidth;
    }

    render(
      g: CanvasRenderingContext2D,
      x: number,
      y: number,
      w: number,
      h: number,
      gridCell: DG.GridCell,
      cellStyle: DG.GridCellStyle,
    ) {
      g.save();
      g.beginPath();
      g.rect(x, y, w, h);
      g.clip();
      g.font = '15px monospace';
      g.textBaseline = 'top';
      const s: string = gridCell.cell.value;
      const textSize = g.measureText(s);
      const subParts = s.split('-');
      x = Math.max(x, x+(w-textSize.width)/2);
      const simplified = !subParts.some((amino, index)=>{

        return amino.length>1&&index!=0&&index!=subParts.length-1;
      });
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

