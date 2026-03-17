import * as DG from 'datagrok-api/dg';
import {RDModule, RDReaction} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {drawErrorCross, drawRdKitReactionToOffscreenCanvas} from '../utils/chem-common-rdkit';
import {USE_RDKIT_REACTION_RENDERER} from '../utils/reactions/consts';

export class RDKitReactionRenderer extends DG.GridCellRenderer {
  rdKitModule: RDModule;
  canvasCounter: number;
  reactionCache: DG.LruCache<String, RDReaction | null> = new DG.LruCache<String, RDReaction | null>();
  reactionRendersCache: DG.LruCache<String, ImageData> = new DG.LruCache<String, ImageData>();
  canvasReused: OffscreenCanvas;

  constructor(rdKitModule: RDModule) {
    super();
    this.rdKitModule = rdKitModule;
    this.canvasCounter = 0;
    this.canvasReused = new OffscreenCanvas(this.defaultWidth, this.defaultHeight);
    this.reactionCache.onItemEvicted = function(obj: any) {
      obj?.delete();
    };
  }

  ensureCanvasSize(w: number, h: number): OffscreenCanvas {
    if (this.canvasReused.width < w || this.canvasReused.height < h)
      this.canvasReused = new OffscreenCanvas(Math.max(this.defaultWidth, w), Math.max(this.defaultHeight, h));
    return this.canvasReused;
  }

  get name(): string {return 'RDKit reaction renderer';}
  get cellType(): string {return 'ChemicalReaction';}
  get defaultWidth() {return 600;}
  get defaultHeight() {return 150;}

  _fetchRxnGetOrCreate(reactionString: string, details: object = {}): RDReaction | null {
    let rxn: RDReaction | null = null;
    try {
      rxn = this.rdKitModule.get_rxn(reactionString, JSON.stringify(details));
    } catch (e) { }
    if (!rxn)
      console.error('Chem | In _fetchMolGetOrCreate: RDKit .get_rxn crashes on a molString: `' + reactionString + '`');
    return rxn;
  }

  _fetchRxn(reactionString: string, details: object = {}): RDReaction | null {
    const name = reactionString + ' || ' + (Object.keys(details).length ? ' || ' + JSON.stringify(details) : '');
    return this.reactionCache.getOrCreate(name, (_: any) => this._fetchRxnGetOrCreate(reactionString, details));
  }

  _rendererGetOrCreate(
    width: number, height: number, reactionString: string): ImageData {
    const rdkitRxn = this._fetchRxn(reactionString);

    const canvas = this.ensureCanvasSize(width, height);//new OffscreenCanvas(width, height);
    const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
    this.canvasCounter++;
    if (rdkitRxn != null)
      drawRdKitReactionToOffscreenCanvas(rdkitRxn, width, height, canvas);

    else {
      // draw a crossed rectangle
      drawErrorCross(ctx, width, height);
    }

    return ctx.getImageData(0, 0, width, height);
  }

  _fetchRender(width: number, height: number, reactionString: string): ImageData {
    const name = width + ' || ' + height + ' || ' + reactionString;
    return this.reactionRendersCache.getOrCreate(name, (_: any) =>
      this._rendererGetOrCreate(width, height, reactionString));
  }

  _drawReaction(x: number, y: number, w: number, h: number, onscreenCanvas: HTMLCanvasElement,
    reactionString: string): void {
    const imageData = this._fetchRender(w, h, reactionString);
        onscreenCanvas.getContext('2d', {willReadFrequently: true})!.putImageData(imageData, x, y);
  }

  /** Render a reaction SMARTS string directly onto an HTML canvas element.
   *  Useful for standalone preview canvases outside of grid cells. */
  renderToCanvas(canvas: HTMLCanvasElement, reactionString: string, width?: number, height?: number): boolean {
    if (!reactionString)
      return false;
    const w = width ?? canvas.width;
    const h = height ?? canvas.height;
    try {
      this._drawReaction(0, 0, w, h, canvas, reactionString);
      return true;
    } catch {
      return false;
    }
  }

  render(g: any, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell): void {
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

// ---- Standalone rendering utility ----

let _sharedRenderer: RDKitReactionRenderer | null = null;

/**
 * Render a reaction SMARTS to an HTML canvas.
 * When USE_RDKIT_REACTION_RENDERER is true, uses the RDKitReactionRenderer class
 * with offscreen canvas + LRU cache. When false, uses simple direct RDKit rendering:
 * get_rxn() → draw_to_canvas() → delete().
 */
export function renderReactionToCanvas(
  rdkit: RDModule, canvas: HTMLCanvasElement, smarts: string, w?: number, h?: number,
): boolean {
  if (!smarts) return false;
  const width = w ?? canvas.width;
  const height = h ?? canvas.height;

  if (USE_RDKIT_REACTION_RENDERER) {
    if (!_sharedRenderer) _sharedRenderer = new RDKitReactionRenderer(rdkit);
    return _sharedRenderer.renderToCanvas(canvas, smarts, width, height);
  }

  // Simple direct RDKit approach
  let rxn: RDReaction | null = null;
  try {
    rxn = rdkit.get_rxn(smarts);
    if (!rxn) return false;
    rxn.draw_to_canvas(canvas, width, height);
    return true;
  } catch {
    return false;
  } finally {
    try {rxn?.delete();} catch {/* noop */}
  }
}
