import * as DG from 'datagrok-api/dg';
import {RDModule, RDReaction} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {drawErrorCross, drawRdKitReactionToOffscreenCanvas} from '../utils/chem-common-rdkit';
import {USE_RDKIT_REACTION_RENDERER} from '../utils/reactions/consts';
import {_convertMolNotation, MALFORMED_MOL_V2000} from '../utils/convert-notation-utils';
import {MESSAGE_MALFORMED} from '../constants';

/** Width in pixels reserved for the arrow drawn between consecutive reaction steps. */
const STEP_ARROW_WIDTH = 30;

/** Maximum number of reaction steps rendered horizontally before wrapping to the next row. */
const MAX_STEPS_PER_ROW = 3;

/**
 * Special branch-boundary delimiter. Inserted by route-formatting code (e.g. the reaction
 * enumerator's `formatRoute`) to mark transitions between independent synthesis branches that
 * later merge — e.g. branch A ends with product P1, branch B starts with fresh reactants —
 * so the renderer doesn't draw a fake "P1 → fresh-reactants" reaction at the seam. The sequence
 * `--**--` is not valid SMILES/SMARTS syntax (`*` is a wildcard atom but must be inside brackets
 * or alone, never as `**`; `--` is not a valid bond expression), so it can't collide with real
 * reaction strings.
 */
export const BRANCH_DELIMITER = '--**--';

/**
 * Split a multi-step reaction string into branches of single-step reactions.
 *
 * Returns a 2-D array: the outer array is one entry per branch (separated by BRANCH_DELIMITER
 * `--**--`); each inner array is the linear sequence of single-step reactions inside that branch.
 *
 * Within a branch, `A.B>>C>>D>>E` becomes `[A.B>>C, C>>D, D>>E]` (consecutive `>>`-separated
 * segments paired into single-step reactions). Single-step strings are returned as a single entry.
 *
 * Examples:
 *   "A>>B"                                    → [["A>>B"]]
 *   "A>>B>>C"                                 → [["A>>B", "B>>C"]]                      (1 branch, linear)
 *   "A>>B--**--C>>D"                          → [["A>>B"], ["C>>D"]]                    (2 branches)
 *   "A>>B--**--C>>D>>E--**--B.E>>F"           → [["A>>B"], ["C>>D","D>>E"], ["B.E>>F"]] (3 branches)
 *
 * The renderer draws step-arrows between consecutive steps within a branch, and a different
 * "branch separator" between branches — so a convergent route doesn't get a fake reaction arrow
 * drawn at the seam where one synthesis path ends and the next begins.
 */
export function parseMultiStepReaction(reactionString: string): string[][] {
  const branchStrings = reactionString.split(BRANCH_DELIMITER);
  const branches: string[][] = [];
  for (const branchStr of branchStrings) {
    const parts = branchStr.split('>>');
    const branchSteps: string[] = [];
    if (parts.length <= 2) {
      const trimmed = branchStr.trim();
      if (trimmed) branchSteps.push(trimmed);
    } else {
      for (let i = 0; i < parts.length - 1; i++) {
        const left = parts[i].trim();
        const right = parts[i + 1].trim();
        if (left && right) branchSteps.push(`${left}>>${right}`);
      }
    }
    if (branchSteps.length > 0) branches.push(branchSteps);
  }
  return branches.length > 0 ? branches : [[reactionString]];
}

/**
 * Flatten the branches structure into a single ordered list of step strings. Whether a gap
 * between two consecutive steps was a within-branch continuation or a branch boundary is no
 * longer carried here — both are rendered with the same step arrow. The branch structure still
 * matters at parse time because it prevents fake "prev_product → fresh_reactants" intermediate
 * steps from being introduced at branch seams; visually the arrows are uniform.
 */
function flattenBranches(branches: string[][]): string[] {
  const steps: string[] = [];
  for (const branch of branches) for (const s of branch) steps.push(s);
  return steps;
}

/** Color used for the step arrows and labels between reaction steps. */
const STEP_ARROW_COLOR = '#2083D5';

/**
 * Draw a right-pointing arrow on the canvas to visually connect consecutive
 * reaction steps. The arrow is drawn in STEP_ARROW_COLOR with a "Step N" label above it.
 */
function drawStepArrow(
  ctx: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
  stepNumber: number,
): void {
  const cx = x + w / 2;
  const cy = y + h / 2;
  const arrowLen = w * 0.6;
  const headSize = Math.min(6, h * 0.15);
  const x1 = cx - arrowLen / 2;
  const x2 = cx + arrowLen / 2;

  ctx.save();
  ctx.strokeStyle = STEP_ARROW_COLOR;
  ctx.fillStyle = STEP_ARROW_COLOR;
  ctx.lineWidth = 1.5;

  // "Step N" label above the arrow
  const fontSize = Math.max(8, Math.min(11, h * 0.1));
  ctx.font = `${fontSize}px sans-serif`;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText(`Step ${stepNumber}`, cx, cy - 6);

  // Shaft
  ctx.beginPath();
  ctx.moveTo(x1, cy);
  ctx.lineTo(x2, cy);
  ctx.stroke();

  // Arrowhead
  ctx.beginPath();
  ctx.moveTo(x2, cy);
  ctx.lineTo(x2 - headSize, cy - headSize);
  ctx.lineTo(x2 - headSize, cy + headSize);
  ctx.closePath();
  ctx.fill();

  ctx.restore();
}

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

    const canvas = this.ensureCanvasSize(width, height);
    const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
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

  private _smilesToSmartsCache: DG.LruCache<string, string> = new DG.LruCache<string, string>(500);

  private _reactionToSmarts(rxnString: string): string {
    try {
      const parts = rxnString.split('>>');
      return parts.map((part) => {
        return this._smilesToSmartsCache.getOrCreate(part, (p) => {
          try {
            if (p.includes('M  END') || p.includes('#') || p.includes('$'))
              return p;
            const res =
              _convertMolNotation(p, DG.chem.Notation.Unknown, DG.chem.Notation.Smarts, this.rdKitModule, false);
            return !!res && res !== MESSAGE_MALFORMED ? res : p;
          } catch (e) {
            return p;
          }
        });
      }).join('>>');
    } catch (e) {
      return rxnString;
    }
  }

  _drawReaction(x: number, y: number, w: number, h: number, onscreenCanvas: HTMLCanvasElement,
    reactionString: string): void {
    // reactions can be in smarts or smiles. rdkit will treat both as smarts and will render smiles in a very bad way.
    const reactionSmarts = this._reactionToSmarts(reactionString);

    const imageData = this._fetchRender(w, h, reactionSmarts);
        onscreenCanvas.getContext('2d', {willReadFrequently: true})!.putImageData(imageData, x, y);
  }

  /**
   * Render a multi-step reaction onto a canvas by splitting it into individual
   * steps laid out horizontally (up to MAX_STEPS_PER_ROW per row, then wrapping).
   * A step arrow is drawn between every pair of consecutive steps. The branch structure from
   * parseMultiStepReaction is preserved at parse time so no fake intermediate steps slip in
   * at branch seams; visually the arrows are uniform.
   */
  _drawMultiStepReaction(
    x: number, y: number, w: number, h: number,
    canvas: HTMLCanvasElement, branches: string[][],
  ): void {
    const steps = flattenBranches(branches);
    const numRows = Math.ceil(steps.length / MAX_STEPS_PER_ROW);
    const rowH = Math.floor(h / numRows);
    const isLastRow = (row: number) => row === numRows - 1;
    const dpr = window.devicePixelRatio;
    const arrowWidth = STEP_ARROW_WIDTH;
    for (let row = 0; row < numRows; row++) {
      const rowStart = row * MAX_STEPS_PER_ROW;
      const rowEnd = Math.min(rowStart + MAX_STEPS_PER_ROW, steps.length);
      const stepsInRow = rowEnd - rowStart;

      // Reserve arrow slots between steps, plus one at the end if the row wraps.
      const hasWrapArrow = !isLastRow(row);
      const arrowCount = stepsInRow - 1 + (hasWrapArrow ? 1 : 0);
      const totalArrowW = arrowCount * arrowWidth * dpr;
      const stepW = Math.floor((w - totalArrowW) / stepsInRow);
      const rowY = y + row * rowH;
      const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
      for (let i = 0; i < stepsInRow; i++) {
        const stepX = x + i * (stepW + arrowWidth * dpr);
        this._drawReaction(stepX, rowY, stepW, rowH, canvas, steps[rowStart + i]);

        // Draw a step arrow after this step (if there's another to follow).
        const globalIdx = rowStart + i;
        if (globalIdx < steps.length - 1) {
          const arrowX = stepX + stepW;
          drawStepArrow(ctx, arrowX / dpr, rowY / dpr, arrowWidth, rowH / dpr, globalIdx + 2);
        }
      }
    }
  }

  /** Render a reaction SMARTS string directly onto an HTML canvas element.
   *  Supports multi-step reactions (multiple ">>" separators) by splitting
   *  into individual steps and rendering them side by side with arrows. */
  renderToCanvas(canvas: HTMLCanvasElement, reactionString: string, width?: number, height?: number): boolean {
    if (!reactionString)
      return false;
    const w = width ?? canvas.width;
    const h = height ?? canvas.height;
    try {
      const branches = parseMultiStepReaction(reactionString);
      const totalSteps = branches.reduce((n, b) => n + b.length, 0);
      if (totalSteps === 1)
        this._drawReaction(0, 0, w, h, canvas, branches[0][0]);
      else
        this._drawMultiStepReaction(0, 0, w, h, canvas, branches);
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
    if (totalSteps === 1) {
      // Single-step reaction — render directly
      this._drawReaction(x, y, w, h, g.canvas, branches[0][0]);
    } else {
      // Multi-step reaction — lay out steps horizontally with arrows / branch separators.
      this._drawMultiStepReaction(x, y, w, h, g.canvas, branches);
    }
  }
}

// ---- Standalone rendering utility ----

let _sharedRenderer: RDKitReactionRenderer | null = null;

/**
 * Render a reaction SMARTS to an HTML canvas.
 * Supports both single-step ("A>>B") and multi-step ("A>>B>>C>>D") reactions.
 *
 * When USE_RDKIT_REACTION_RENDERER is true, uses the RDKitReactionRenderer class
 * with offscreen canvas + LRU cache. When false, uses simple direct RDKit rendering.
 *
 * For multi-step reactions in direct mode, each step is rendered separately
 * side by side with arrows connecting them.
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

  // Direct RDKit rendering — handle multi-step by splitting
  const branches = parseMultiStepReaction(smarts);
  const steps = flattenBranches(branches);

  if (steps.length === 1) {
    // Single step — simple direct rendering
    let rxn: RDReaction | null = null;
    try {
      rxn = rdkit.get_rxn(steps[0]);
      if (!rxn) return false;
      rxn.draw_to_canvas(canvas, width, height);
      return true;
    } catch {
      return false;
    } finally {
      try {rxn?.delete();} catch {/* noop */}
    }
  }

  // Multi-step — render each step into its own region of the canvas
  const ctx = canvas.getContext('2d')!;
  const arrowCount = steps.length - 1;
  const dpr = window.devicePixelRatio || 1;
  const arrowWidth = STEP_ARROW_WIDTH * dpr;
  const totalArrowW = arrowCount * arrowWidth;
  const stepW = Math.floor((width - totalArrowW) / steps.length);

  for (let i = 0; i < steps.length; i++) {
    let rxn: RDReaction | null = null;
    try {
      rxn = rdkit.get_rxn(steps[i]);
      if (!rxn) continue;

      // Render this step into a temporary canvas, then composite onto the target
      const tmpCanvas = document.createElement('canvas');
      tmpCanvas.width = stepW;
      tmpCanvas.height = height;
      rxn.draw_to_canvas(tmpCanvas, stepW, height);

      const stepX = i * (stepW + arrowWidth);
      ctx.drawImage(tmpCanvas, stepX, 0);

      // Draw a step arrow after each step except the last.
      if (i < steps.length - 1)
        drawStepArrow(ctx, stepX + stepW, 0, arrowWidth, height, i + 2);
    } catch {
      // skip failed step
    } finally {
      try {rxn?.delete();} catch {/* noop */}
    }
  }
  return true;
}
