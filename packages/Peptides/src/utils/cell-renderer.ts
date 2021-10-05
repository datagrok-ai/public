import {ChemPalette} from './chem-palette';
import * as DG from 'datagrok-api/dg';

const cp = ChemPalette.getDatagrok();

function getColorPivot(c = ''): [string, number] {
  if (c.length == 1 || c.at(1) == '(') {
    const amino = c.at(0)?.toUpperCase()!;
    return [amino in cp ? cp[amino] : 'rgb(77,77,77)', 1];
  }
  if (c.at(0) == 'd' && c.at(1)! in cp) {
    if (c.length == 2 || c.at(2) == '(') {
      const amino = c.at(1)?.toUpperCase()!;
      return [amino in cp ? cp[amino] : 'rgb(77,77,77)', 2];
    }
  }
  if (c.substr(0, 3) in ChemPalette.AAFullNames) {
    if (c.length == 3 || c.at(3) == '(') {
      const amino = ChemPalette.AAFullNames[c.substr(0, 3)];
      return [amino in cp ? cp[amino] : 'rgb(77,77,77)', 3];
    }
  }
  if (c.at(0)?.toLowerCase() == c.at(0)) {
    if (c.substr(1, 3) in ChemPalette.AAFullNames) {
      if (c.length == 4 || c.at(4) == '(') {
        const amino = ChemPalette.AAFullNames[c.substr(1, 3)];
        return [amino in cp ? cp[amino] : 'rgb(77,77,77)', 4];
      }
    }
  }
  return ['rgb(77,77,77)', 0];
  //return c ? DG.Color.toRgb(this.colorScale(c)) : 'rgb(127,127,127)'
};

function printLeftCentered(
  x: number,
  y: number,
  w: number,
  h: number,
  g: CanvasRenderingContext2D,
  s: string,
  color = 'rgb(77,77,77)',
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
    g.fillStyle = 'rgb(77,77,77)';
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
    g.fillStyle = 'rgb(77,77,77)';
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
      const [color, pivot] = getColorPivot(s);
      printLeftCentered(x, y, w, h, g, s, color, pivot, false, true);
      g.restore();
    }
}

export class AlignedSequenceCellRenderer extends DG.GridCellRenderer {
    private maxCellWidth = 270;

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
      const subParts = s.split('-');
      subParts.forEach((amino: string, index) => {
        const [color, pivot] = getColorPivot(amino);
        g.fillStyle = 'rgb(77,77,77)';
        if (index + 1 < subParts.length) {
          amino += '-';
        }
        x = printLeftCentered(x, y, w, h, g, amino, color, pivot, true);
      });
      g.restore();
    }
}

