import {ChemPalette} from './chem-palette';
import * as DG from 'datagrok-api/dg';


export class AminoAcidsCellRenderer extends DG.GridCellRenderer {
    private maxCellWidth = 11;
    private fontSize = 15;
    private fontSizeSide = 35;
    private spacing = 0.8;

    get name() {
      return 'aminoAcidsCR';
    }

    get cellType() {
      return 'aminoAcids';
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
      g.font = 'bold 15px monospace';
      g.textBaseline = 'top';
      const s = gridCell.cell.value ? gridCell.cell.value : '-';
      const cp = ChemPalette.get_datagrok();

      cp['-'] = `rgb(77, 77, 77)`;

      if (s in cp) {
        g.font = `bold ${this.fontSize}px monospace`;
        g.fillStyle = cp[s];
        const textSize = g.measureText(s);
        g.fillText(s, x + (w - textSize.width)/2, y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent)/2);
        return;
      }
      g.font = `${this.fontSizeSide / s.length}px monospace`;
      for (let i = 0; i < s.length; i++) {
        g.fillStyle = `rgb(77, 77, 77)`;
        g.fillText(s[i], x + 3, y + i * this.fontSizeSide * this.spacing / s.length);
      }
    }
}

export class AlignedSequenceCellRenderer extends DG.GridCellRenderer {
    private maxCellWidth = 270;
    private spacing = 9;
    private fontSize = 18;

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
      g.font = 'bold 15px monospace';
      g.textBaseline = 'top';
      const s = gridCell.cell.value;

      const cp = ChemPalette.get_datagrok();

      let xPos = 0;
      let yPos = 0;
      let skipped = 0;
      let odd = true;

      g.font = `${this.fontSize}}px monospace`;
      for (let i = 0; i < Math.min(s.length, 3); i++) {
        xPos = 2 + x + i * this.spacing;
        yPos = y + 6;
        g.fillStyle = `rgb(128, 128, 128)`;
        g.fillText(s[i], xPos, yPos);
      }
      g.font = `bold ${this.fontSize}}px monospace`;
      for (let i = 3; i < s.length - 4; i++) {
        xPos = 3 + x + (i - skipped) * this.spacing;
        yPos = y + 6;
        if (s.charAt(i) in cp) {
          g.fillStyle = cp[s.charAt(i)];
          g.fillText(s[i], xPos, yPos);
          odd = true;
          continue;
        }
        if (odd) {
          skipped += 1;
          odd = false;
          continue;
        }
        g.fillStyle = `rgb(128, 128, 128)`;
        g.fillText(s[i], xPos, yPos);
      }
      g.font = `${this.fontSize}}px monospace`;
      for (let i = Math.max(0, s.length - 4); i < s.length; i++) {
        xPos = 4 + x + (i - skipped) * this.spacing;
        yPos = y + 6;
        g.fillStyle = `rgb(128, 128, 128)`;
        g.fillText(s[i], xPos, yPos);
      }
    }
}

