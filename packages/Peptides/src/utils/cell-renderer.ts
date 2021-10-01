import {ChemPalette} from './chem-palette';
import * as DG from 'datagrok-api/dg';

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
      const s = gridCell.cell.value ? gridCell.cell.value : '-';
      const cp = ChemPalette.getDatagrok();

      cp['-'] = `rgb(77, 77, 77)`;
      g.fillStyle =`rgb(77, 77, 77)`;
      if (s in cp) {
        g.fillStyle = cp[s];
      }
      const textSize = g.measureText(s);
      g.fillText(
        s,
        x + (w - textSize.width) / 2,
        y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2,
      );
      g.restore();
    }
}

export class AlignedSequenceCellRenderer extends DG.GridCellRenderer {
    private maxCellWidth = 270;
    private spacing = 9;
    private fontSize = 18;
    private cp = ChemPalette.getDatagrok();
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
      g.scale(1, 1);
      g.clearRect(0, 0, h, w);
      g.font = '15px monospace';
      g.textBaseline = 'top';
      const s = gridCell.cell.value;


      let xPos = 0;
      let yPos = 0;
      let skipped = 0;
      let odd = true;

      g.font = `${this.fontSize}}px monospace`;
      for (let i = 0; i < Math.min(s.length, 3); i++) {
        xPos = x+2 + i * this.spacing;
        yPos =y + 6;
        g.fillStyle = `rgb(128, 128, 128)`;
        g.fillText(s[i], xPos, yPos);
      }
      g.font = `${this.fontSize}}px monospace`;
      for (let i = 3; i < s.length - 4; i++) {
        xPos =x+ 3 + (i - skipped) * this.spacing;
        yPos = y+ 6;
        if (s.charAt(i) in this.cp) {
          g.fillStyle = this.cp[s.charAt(i)];
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
        xPos = x+ 4 + (i - skipped) * this.spacing;
        yPos = y+ 6;
        g.fillStyle = `rgb(128, 128, 128)`;
        g.fillText(s[i], xPos, yPos);
      }
      g.restore();
    }
}

