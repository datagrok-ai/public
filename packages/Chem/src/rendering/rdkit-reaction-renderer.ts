import * as DG from 'datagrok-api/dg';
import {RDModule, RDReaction} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {drawErrorCross, drawRdKitReactionToOffscreenCanvas} from '../utils/chem-common-rdkit';
import {USE_RDKIT_REACTION_RENDERER} from '../utils/reactions/consts';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {MESSAGE_MALFORMED} from '../constants';
import {prepareReactionDepiction} from './reaction-depiction';

const MAX_STEPS_PER_ROW = 3;
const MIN_STEP_WIDTH = 240;

/** Separates independent steps/branches without inventing a reaction at their boundary. */
export const BRANCH_DELIMITER = '--**--';

/** A.B>>C>>D becomes [[A.B>>C, C>>D]]; --**-- starts a separate branch. */
export function parseMultiStepReaction(reactionString: string): string[][] {
  const branchStrings = reactionString.split(BRANCH_DELIMITER);
  const branches: string[][] = [];
  for (const branchStr of branchStrings) {
    const parts = branchStr.split('>>');
    const branchSteps: string[] = [];
    if (parts.length <= 2) {
      const trimmed = branchStr.includes('M  END') ? branchStr : branchStr.trim();
      if (trimmed) branchSteps.push(trimmed);
    } else {
      for (let i = 0; i < parts.length - 1; i++) {
        const left = parts[i].includes('M  END') ? parts[i] : parts[i].trim();
        const right = parts[i + 1].includes('M  END') ? parts[i + 1] : parts[i + 1].trim();
        if (left && right) branchSteps.push(`${left}>>${right}`);
      }
    }
    if (branchSteps.length > 0) branches.push(branchSteps);
  }
  return branches.length > 0 ? branches : [[reactionString]];
}

export class RDKitReactionRenderer extends DG.GridCellRenderer {
  rdKitModule: RDModule;
  canvasCounter: number;
  reactionCache: DG.LruCache<string, RDReaction | null> = new DG.LruCache<string, RDReaction | null>();
  reactionRendersCache: DG.LruCache<string, ImageData> = new DG.LruCache<string, ImageData>();
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
    } catch { }
    if (!rxn)
      console.error('Chem | RDKit could not parse reaction: `' + reactionString + '`');
    return rxn;
  }

  _fetchRxn(reactionString: string, details: object = {}): RDReaction | null {
    const name = reactionString + ' || ' + (Object.keys(details).length ? ' || ' + JSON.stringify(details) : '');
    return this.reactionCache.getOrCreate(name, (_: any) => this._fetchRxnGetOrCreate(reactionString, details));
  }

  _rendererGetOrCreate(
    width: number, height: number, reactionString: string): ImageData {
    const rdkitRxn = this._fetchRxn(reactionString);

    const canvas = this.ensureCanvasSize(width, height);
    const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
    ctx.clearRect(0, 0, width, height);
    this.canvasCounter++;
    if (rdkitRxn != null)
      drawRdKitReactionToOffscreenCanvas(rdkitRxn, width, height, canvas);
    else
      drawErrorCross(ctx, width, height);

    return ctx.getImageData(0, 0, width, height);
  }

  _fetchRender(width: number, height: number, reactionString: string): ImageData {
    const name = width + ' || ' + height + ' || ' + reactionString;
    return this.reactionRendersCache.getOrCreate(name, (_: any) =>
      this._rendererGetOrCreate(width, height, reactionString));
  }

  private _depictionsCache = new DG.LruCache<string, {steps: string[], molblocks: Map<string, string>}>(500);

  private _prepareSteps(steps: string[]): string[] {
    return this._depictionsCache.getOrCreate(steps.join(BRANCH_DELIMITER), () => {
      const molblocks = new Map<string, string>();
      return {
        steps: steps.map((s) => prepareReactionDepiction(this.rdKitModule, s, molblocks) ?? this._reactionToSmarts(s)),
        molblocks,
      };
    }).steps;
  }

  private _reactionToSmarts(rxnString: string): string {
    if (rxnString.includes('|') || rxnString.includes('M  END'))
      return rxnString;
    return rxnString.split(/(?<!-)>/).map((part) => {
      if (!part || part.includes('#') || part.includes('$'))
        return part;
      const res = _convertMolNotation(part, DG.chem.Notation.Unknown, DG.chem.Notation.Smarts, this.rdKitModule, false);
      return res && res !== MESSAGE_MALFORMED ? res : part;
    }).join('>');
  }

  _drawReaction(x: number, y: number, w: number, h: number, onscreenCanvas: HTMLCanvasElement,
    reactionString: string): void {
    if (w < 1 || h < 1)
      return;
    const imageData = this._fetchRender(Math.floor(w), Math.floor(h), reactionString);
    onscreenCanvas.getContext('2d', {willReadFrequently: true})!.putImageData(imageData, Math.round(x), Math.round(y));
  }

  _drawMultiStepReaction(
    x: number, y: number, w: number, h: number,
    canvas: HTMLCanvasElement, branches: string[][], pixelRatio = window.devicePixelRatio,
  ): void {
    const steps = this._prepareSteps(branches.flat());
    const columns = Math.min(MAX_STEPS_PER_ROW, Math.max(1, Math.floor(w / pixelRatio / MIN_STEP_WIDTH)));
    const numRows = Math.ceil(steps.length / columns);
    const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
    const style = getComputedStyle(canvas);
    ctx.save();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.beginPath();
    ctx.rect(x, y, w, h);
    ctx.clip();
    const labelFont = `${11 * pixelRatio}px sans-serif`;
    ctx.font = labelFont;
    const prefixWidth = ctx.measureText('Step ').width;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    ctx.lineWidth = pixelRatio;
    try {
      for (let row = 0; row < numRows; row++) {
        const rowStart = row * columns;
        const stepsInRow = Math.min(columns, steps.length - rowStart);
        const rowY = y + Math.floor(row * h / numRows);
        const rowH = Math.floor((row + 1) * h / numRows) - Math.floor(row * h / numRows);
        const headerH = rowH >= 48 * pixelRatio ? 20 * pixelRatio : 0;
        const gap = Math.min(12 * pixelRatio, w / stepsInRow / 10);
        for (let i = 0; i < stepsInRow; i++) {
          const stepX = x + Math.floor(i * w / stepsInRow);
          const stepW = Math.floor((i + 1) * w / stepsInRow) - Math.floor(i * w / stepsInRow);
          this._drawReaction(stepX + gap / 2, rowY + headerH, stepW - gap, rowH - headerH,
            canvas, steps[rowStart + i]);
          if (headerH) {
            ctx.fillStyle = style.getPropertyValue('--blue-1').trim() || style.color;
            ctx.fillText('Step ', stepX + gap, rowY + 4 * pixelRatio);
            ctx.font = `bold ${labelFont}`;
            ctx.fillText(`${rowStart + i + 1}`, stepX + gap + prefixWidth, rowY + 4 * pixelRatio);
            ctx.font = labelFont;
          }
          ctx.strokeStyle = style.getPropertyValue('--grey-2').trim() || style.color;
          ctx.beginPath();
          if (i > 0) {
            ctx.moveTo(stepX, rowY + 4 * pixelRatio);
            ctx.lineTo(stepX, rowY + rowH - 4 * pixelRatio);
          }
          if (row > 0) {
            ctx.moveTo(stepX + gap, rowY);
            ctx.lineTo(stepX + stepW - gap, rowY);
          }
          ctx.stroke();
        }
      }
    } finally {
      ctx.restore();
    }
  }

  /** Render a reaction SMARTS string directly onto an HTML canvas element.
   *  Supports multi-step reactions (multiple ">>" separators) by splitting
   *  into numbered panels that retain only the chemical reaction arrows. */
  renderToCanvas(canvas: HTMLCanvasElement, reactionString: string, width?: number, height?: number): boolean {
    if (!reactionString)
      return false;
    const w = width ?? canvas.width;
    const h = height ?? canvas.height;
    if (w < 1 || h < 1)
      return false;
    try {
      const ctx = canvas.getContext('2d')!;
      ctx.save();
      ctx.setTransform(1, 0, 0, 1, 0, 0);
      ctx.clearRect(0, 0, w, h);
      ctx.restore();
      const branches = parseMultiStepReaction(reactionString);
      const totalSteps = branches.reduce((n, b) => n + b.length, 0);
      if (totalSteps === 1)
        this._drawReaction(0, 0, w, h, canvas, this._prepareSteps([branches[0][0]])[0]);
      else {
        this._drawMultiStepReaction(0, 0, w, h, canvas, branches,
          canvas.clientWidth > 0 ? canvas.width / canvas.clientWidth : 1);
      }
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

    const branches = parseMultiStepReaction(reactionString);
    const totalSteps = branches.reduce((n, b) => n + b.length, 0);
    if (totalSteps === 1)
      this._drawReaction(x, y, w, h, g.canvas, this._prepareSteps([branches[0][0]])[0]);
    else
      this._drawMultiStepReaction(x, y, w, h, g.canvas, branches);
  }
}

let _sharedRenderer: RDKitReactionRenderer | null = null;

/** Draw a single reaction or a route in the same layout used by the grid. */
export function renderReactionToCanvas(
  rdkit: RDModule, canvas: HTMLCanvasElement, smarts: string, w?: number, h?: number,
): boolean {
  if (!smarts) return false;
  const width = w ?? canvas.width;
  const height = h ?? canvas.height;
  const steps = parseMultiStepReaction(smarts).flat();

  if (USE_RDKIT_REACTION_RENDERER || steps.length > 1) {
    if (!_sharedRenderer || _sharedRenderer.rdKitModule !== rdkit)
      _sharedRenderer = new RDKitReactionRenderer(rdkit);
    return _sharedRenderer.renderToCanvas(canvas, smarts, width, height);
  }

  let rxn: RDReaction | null = null;
  try {
    rxn = rdkit.get_rxn(steps[0]);
    if (!rxn) return false;
    rxn.draw_to_canvas(canvas, width, height);
    return true;
  } catch {
    return false;
  } finally {
    rxn?.delete();
  }
}
