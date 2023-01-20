import * as DG from 'datagrok-api/dg';
import { RDModule, Reaction } from "@datagrok-libraries/chem-meta/src/rdkit-api";
import { drawErrorCross, drawRdKitReactionToOffscreenCanvas } from '../utils/chem-common-rdkit';

export class RDKitReactionRenderer extends DG.GridCellRenderer {

    rdKitModule: RDModule;
    canvasCounter: number;
    reactionCache: DG.LruCache<String, Reaction | null> = new DG.LruCache<String, Reaction | null>();
    reactionRendersCache: DG.LruCache<String, ImageData> = new DG.LruCache<String, ImageData>();
    canvasReused: OffscreenCanvas;

    constructor(rdKitModule: RDModule) {
        super();
        this.rdKitModule = rdKitModule;
        this.canvasCounter = 0;
        this.canvasReused = new OffscreenCanvas(this.defaultWidth, this.defaultHeight);
        this.reactionCache.onItemEvicted = function (obj: any) {
            obj?.delete();
        };
    }

    ensureCanvasSize(w: number, h: number): OffscreenCanvas {
        if (this.canvasReused.width < w || this.canvasReused.height < h)
            this.canvasReused = new OffscreenCanvas(Math.max(this.defaultWidth, w), Math.max(this.defaultHeight, h));
        return this.canvasReused;
    }

    get name(): string { return 'RDKit reaction renderer'; }
    get cellType(): string { return 'ChemicalReaction'; }
    get defaultWidth() { return 600; }
    get defaultHeight() { return 150; }

    _fetchRxnGetOrCreate(reactionString: string, details: object = {}): Reaction | null {
        let rxn: Reaction | null = null;
        try {
            rxn = this.rdKitModule.get_rxn(reactionString, JSON.stringify(details));
        } catch (e) { }
        if (!rxn)
            console.error('Chem | In _fetchMolGetOrCreate: RDKit .get_rxn crashes on a molString: `' + reactionString + '`');
        return rxn;
    }

    _fetchRxn(reactionString: string, details: object = {}): Reaction | null {
        const name = reactionString + ' || ' + (Object.keys(details).length ? ' || ' + JSON.stringify(details) : '');
        return this.reactionCache.getOrCreate(name, (_: any) => this._fetchRxnGetOrCreate(reactionString, details));
    }

    _rendererGetOrCreate(
        width: number, height: number, reactionString: string): ImageData {
        const rdkitRxn = this._fetchRxn(reactionString, { useSmiles: true });

        const canvas = this.ensureCanvasSize(width, height);//new OffscreenCanvas(width, height);
        const ctx = canvas.getContext('2d', { willReadFrequently: true })!;
        this.canvasCounter++;
        if (rdkitRxn != null) {
            drawRdKitReactionToOffscreenCanvas(rdkitRxn, width, height, canvas);
        }
        else {
            // draw a crossed rectangle
            drawErrorCross(ctx, width, height);
        }

        return ctx.getImageData(0, 0, width, height);
    }

    _fetchRender(width: number, height: number, reactionString: string): ImageData {
        const name = width + ' || ' + height + ' || ' + reactionString;
        return this.reactionRendersCache.getOrCreate(name, (_: any) => this._rendererGetOrCreate(width, height, reactionString));
    }

    _drawReaction(x: number, y: number, w: number, h: number, onscreenCanvas: HTMLCanvasElement, reactionString: string): void {
        const imageData = this._fetchRender(w, h, reactionString);
        onscreenCanvas.getContext('2d', { willReadFrequently: true })!.putImageData(imageData, x, y);
    }

    render(g: any, x: number, y: number, w: number, h: number,
        gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
        const reactionString = gridCell.cell.value;
        if (reactionString == null || reactionString === '')
            return;

        const r = window.devicePixelRatio;
        x = r * x; y = r * y;
        w = r * w; h = r * h;

        const reactionParts = reactionString.length > 130 ? reactionString.split('>>') : [reactionString];
        if (reactionParts.length > 1) {
            reactionParts[0] = `${reactionParts[0]}>>`;
            reactionParts[1] = `>>${reactionParts[1]}`;
            h = h/2;
        }
        let yUpdated = y;
        reactionParts.forEach((part: string) => {
            this._drawReaction(x, yUpdated, w, h, g.canvas, part);
            yUpdated = yUpdated + h;
        });
    }


}
  