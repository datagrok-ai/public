import {ChemPalette} from './chem-palette';
import * as DG from 'datagrok-api/dg';

const cp = ChemPalette.getDatagrok();
function getColor(c = '') {
  if (c.length == 1 || c.at(1) == '(' || c.at(0)?.toLowerCase() == 'd') {
    const amino = c.at(0)?.toUpperCase()!;
    return amino in cp ? cp[amino] : 'rgb(77,77,77)';
  }
  return 'rgb(77,77,77)';
  //return c ? DG.Color.toRgb(this.colorScale(c)) : 'rgb(127,127,127)'
};

function printLeftCentered(
  x:number,
  y:number,
  w:number,
  h:number,
  g :CanvasRenderingContext2D,
  s:string,
  color = 'rgb(77,77,77)',
  pivot:number = 0,
  left = false,
) {
  const textSize = g.measureText(s);
  const indent = 5;
  const colorPart = pivot == -1 ? s.substring(0) : s.substring(0, pivot);
  const grayPart = pivot == -1 ? '' : s.substr(pivot);

  const colorTextSize = g.measureText(colorPart);
  if (left || textSize.width > w) {
    g.fillStyle = color;
    g.fillText(
      colorPart,
      x + indent,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    g.fillStyle = getColor();
    g.fillText(
      grayPart,
      x + indent + colorTextSize.width,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    return x + colorTextSize.width + g.measureText(grayPart).width; ;
  } else {
    g.fillStyle = color;
    g.fillText(
      colorPart,
      x + (w - textSize.width) / 2,
      y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
    );
    g.fillStyle = getColor();
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
      const color = getColor(s);
      if (s.at(0)?.toLowerCase() == 'd') {
        const modInd = s.indexOf('(');
        printLeftCentered(x, y, w, h, g, s, color, modInd);
      } else if (s.at(1) == '(') {
        printLeftCentered(x, y, w, h, g, s, color, 1);
      } else {
        printLeftCentered(x, y, w, h, g, s, color, -1);
      }
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
      const s:string = gridCell.cell.value;
      const subParts = s.split('-');
      subParts.forEach((amino:string, index) =>{
        const color = getColor(amino);
        g.fillStyle = 'rgb(77,77,77)';
        if (index+1 < subParts.length) {
          amino += '-';
        }
        if (amino.at(0)?.toLowerCase() == 'd') {
          const modInd = amino.indexOf('(');
          x = printLeftCentered(x, y, w, h, g, amino, color, modInd, true);
        } else if (amino.at(1) == '(') {
          x =printLeftCentered(x, y, w, h, g, amino, color, 1, true);
        } else {
          x =printLeftCentered(x, y, w, h, g, amino, color, -1, true);
        }
      });
      g.restore();
    }
}

