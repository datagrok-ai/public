import {ChemPalette} from './chem-palette';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class AlignedSequenceCellRenderer extends DG.GridCellRenderer {
    get name() {
        return 'alignedSequenceCR';
    }

    get cellType() {
        return 'alignedSequence';
    }

    render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
        g.font = 'bold 15px monospace';
        g.textBaseline = 'top';
        let s = gridCell.cell.value;

        let cp = ChemPalette.get_datagrok()
        g.font = '13px monospace';
        for (let i = 0; i < Math.max(s.length, 3); i++) {
            let xPos = 2 + x + i * 7
            let yPos = y + 6
            g.fillStyle = `rgb(128, 128, 128)`;
            g.fillText(s[i], xPos, yPos);
        }
        for (let i = 3; i < s.length - 4; i++) {
            let xPos = 2 + x + i * 7
            let yPos = y + 6
            if (s.charAt(i) in cp) {
                g.fillStyle = cp[s.charAt(i)];
                g.font = '18px monospace';
                g.fillText('█', xPos, yPos);
                g.fillStyle = 'black';
                g.strokeStyle = 'white';
                g.lineWidth = 0.8;
                g.font = 'bold 15px monospace';
                g.strokeText(s[i], xPos+0.5, yPos);

                g.fillText(s[i], xPos+0.5, yPos);
            } else {
                g.fillStyle = `rgb(128, 128, 128)`;
                g.fillText(s[i], xPos, yPos);
            }
        }
        g.font = '13px monospace';
        for (let i = Math.max(0, s.length - 4); i < s.length; i++) {
            let xPos = 2 + x + i * 7
            let yPos = y + 6
            g.fillStyle = `rgb(128, 128, 128)`;
            g.fillText(s[i], xPos, yPos);
        }
    }
}