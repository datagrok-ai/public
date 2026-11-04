import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {monomerToShort} from './macromolecule/utils';
import {ILogger} from './logger';

const undefinedColor = 'rgb(100,100,100)';
const grayColor = '#808080';
const blackColor = 'rgb(0,0,0)';
const monomerToShortFunction: (amino: string, maxLengthOfMonomer: number) => string = monomerToShort;

export enum TAGS {
  maxMonomerLength = '.mm.cellRenderer.maxMonomerLength',
}

export enum MONOMER_RENDERER_TAGS {
  applyToBackground = '.m.cellRenderer.applyToBackground',
}

export enum DrawStyle {
  MSA = 'MSA',
  classic = 'classic',
}

export const PrintOptionsDefaults = new class {
  /** Color of text to print */ color: string = undefinedColor;
  /** */ pivot: number = 0;
  /** Is left aligned */ left: boolean = false;
  /** Transparency rate where 1.0 is fully transparent */ transparencyRate: number = 0;
  /** Monomer's separator, if specified */ separator: string = '';
  /** Is checker if element last or not */ last: boolean = false;
  /** MSA - for aligned, classic - for other seq */ drawStyle: DrawStyle = DrawStyle.classic;
  /** Max word lengths per position (index of list) */ maxWord: number[] = [];
  /** Word index we currently draw */ wordIdx: number = 0;
  /**  */ gridCell: DG.GridCell | null = null;
  /** Reference sequence for diff mode */ referenceSequence: string[] | null = null;
  /** Max length of a monomer */ maxLengthOfMonomer: number | null = null;
  /** Map of monomers' text sizes */ monomerTextSizeMap: { [key: string]: TextMetrics } = {};
  /** */ logger?: ILogger = undefined;
  selectedPosition?: number = undefined;
  // Multiline stuff
  isMultiLineContext: boolean = false;
  lineNumber: number = 0;
  /** {@link gridCell}'s column temps as {@link getPrintColumnTemps} reads them, computed once per cell by a caller
   * that prints many monomers of the cell; read from {@link gridCell} per call when not given */
  columnTemps?: PrintColumnTemps = undefined;
}();

export type PrintOptions = typeof PrintOptionsDefaults;

/** The column temps {@link printLeftOrCentered} draws by */
export type PrintColumnTemps = { colorCode: boolean, compareWithCurrent: boolean, highlightDifference: string };

/** The `color-code`, `compare-with-current` and `highlight-difference` temps of the grid cell's column, with their
 * defaults (also for no cell or no column), as {@link printLeftOrCentered} reads them. */
export function getPrintColumnTemps(gridCell: DG.GridCell | null): PrintColumnTemps {
  let colorCode = true;
  let compareWithCurrent = true;
  let highlightDifference = 'difference';
  if ((gridCell != null) && (gridCell.cell.column != null)) {
    colorCode = gridCell.cell.column.temp['color-code'] ?? true;
    compareWithCurrent = gridCell.cell.column.temp['compare-with-current'] ?? true;
    highlightDifference = gridCell.cell.column.temp['highlight-difference'] ?? 'difference';
  }
  return {colorCode, compareWithCurrent, highlightDifference};
}

/** {...PrintOptionsDefaults, ...options} without the generic object spread per printed monomer: a given key wins,
 * even with the value undefined, as with the spread. */
function withPrintDefaults(o: Partial<PrintOptions>): PrintOptions {
  const d = PrintOptionsDefaults;
  return {
    color: 'color' in o ? o.color! : d.color,
    pivot: 'pivot' in o ? o.pivot! : d.pivot,
    left: 'left' in o ? o.left! : d.left,
    transparencyRate: 'transparencyRate' in o ? o.transparencyRate! : d.transparencyRate,
    separator: 'separator' in o ? o.separator! : d.separator,
    last: 'last' in o ? o.last! : d.last,
    drawStyle: 'drawStyle' in o ? o.drawStyle! : d.drawStyle,
    maxWord: 'maxWord' in o ? o.maxWord! : d.maxWord,
    wordIdx: 'wordIdx' in o ? o.wordIdx! : d.wordIdx,
    gridCell: 'gridCell' in o ? o.gridCell! : d.gridCell,
    referenceSequence: 'referenceSequence' in o ? o.referenceSequence! : d.referenceSequence,
    maxLengthOfMonomer: 'maxLengthOfMonomer' in o ? o.maxLengthOfMonomer! : d.maxLengthOfMonomer,
    monomerTextSizeMap: 'monomerTextSizeMap' in o ? o.monomerTextSizeMap! : d.monomerTextSizeMap,
    logger: 'logger' in o ? o.logger : d.logger,
    selectedPosition: 'selectedPosition' in o ? o.selectedPosition : d.selectedPosition,
    isMultiLineContext: 'isMultiLineContext' in o ? o.isMultiLineContext! : d.isMultiLineContext,
    lineNumber: 'lineNumber' in o ? o.lineNumber! : d.lineNumber,
    columnTemps: 'columnTemps' in o ? o.columnTemps : d.columnTemps,
  };
}

/** A function that prints a string aligned to left or centered.
 * @param {CanvasRenderingContext2D} g Canvas rendering context
 * @param {string} s Text to print
 * @param {number} x x coordinate
 * @param {number} y y coordinate
 * @param {number} w Width
 * @param {number} h Height
 * @param {PrintOptions} options
 * @return {number} x coordinate to start printing at.*/
// Updated printLeftOrCentered function to use font-based spacing
// Replace the multiline section in your printLeftOrCentered function

export function printLeftOrCentered(g: CanvasRenderingContext2D,
  s: string, x: number, y: number, w: number, h: number, options: Partial<PrintOptions>
): number {
  const opts = withPrintDefaults(options);

  const logPrefix: string = `Bio: printLeftOrCentered()`;

  if (opts.isMultiLineContext) {
    g.textBaseline = 'middle';
    g.textAlign = 'center'; // Center text horizontally

    let colorPart = s;
    if (opts.maxLengthOfMonomer != null)
      colorPart = monomerToShortFunction(colorPart, opts.maxLengthOfMonomer);

    // Apply transparency
    const alpha = Math.max(0.1, 1.0 - (opts.transparencyRate ?? 0.0));
    g.globalAlpha = alpha;

    const metrics = g.measureText(colorPart);
    const textHeight = metrics.fontBoundingBoxAscent + metrics.fontBoundingBoxDescent;
    const centerY = y + (h - textHeight) / 2 + metrics.fontBoundingBoxAscent;

    let fillColor = opts.color ?? undefinedColor;
    if (!fillColor || fillColor === undefinedColor)
      fillColor = blackColor;

    g.fillStyle = fillColor;

    if (opts.selectedPosition === opts.wordIdx + 1) {
      g.save();
      g.fillStyle = 'rgba(60, 177, 115, 0.2)';
      g.fillRect(x, y, w, h); // Fill the entire column bounds
      g.restore();
      g.fillStyle = fillColor;
    }

    // To center, provide the midpoint of the bounds
    g.fillText(colorPart, x + w / 2, centerY);

    // Reset alpha and other properties
    g.globalAlpha = 1.0;
    g.textBaseline = 'top';
    g.textAlign = 'start'; // Reset to default

    return x + w;
  }

  // ... rest of the existing single-line logic remains unchanged
  options.logger?.debug(`${logPrefix}, start`);
  g.textAlign = 'start';
  let colorPart = s.substring(0);
  let grayPart = opts.last ? '' : opts.separator;
  if (opts.drawStyle === DrawStyle.MSA) grayPart = '';

  const {colorCode, compareWithCurrent, highlightDifference} = opts.columnTemps ?? getPrintColumnTemps(opts.gridCell);

  if (opts.referenceSequence) {
    const currentMonomerCanonical = opts.referenceSequence[opts.wordIdx];
    if (compareWithCurrent && (opts.referenceSequence.length > 0) && (highlightDifference === 'difference'))
      opts.transparencyRate = (colorPart == currentMonomerCanonical) ? 0.7 : opts.transparencyRate;
    if (compareWithCurrent && (opts.referenceSequence.length > 0) && (highlightDifference === 'equal'))
      opts.transparencyRate = (colorPart != currentMonomerCanonical) ? 0.7 : opts.transparencyRate;
  }
  if (opts.maxLengthOfMonomer != null)
    colorPart = monomerToShortFunction(colorPart, opts.maxLengthOfMonomer);

  const fullText = colorPart + grayPart;
  opts.monomerTextSizeMap[fullText] ??= g.measureText(fullText);
  let textSize: any = opts.monomerTextSizeMap[fullText];
  opts.monomerTextSizeMap[colorPart] ??= g.measureText(colorPart);
  let maxColorTextSize = opts.monomerTextSizeMap[colorPart].width;
  const dy = h / 2 - (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2 + 1;

  opts.monomerTextSizeMap[grayPart] ??= g.measureText(grayPart);
  const grayPartSize = opts.monomerTextSizeMap[grayPart].width;
  textSize = textSize.width;
  if (opts.drawStyle === DrawStyle.MSA) {
    maxColorTextSize = opts.maxWord[opts.wordIdx];
    textSize = opts.maxWord[opts.wordIdx];
  }

  function draw(dx1: number, dx2: number): void {
    let drawColor = colorCode ? opts.color : blackColor;
    if (opts.selectedPosition === opts.wordIdx + 1) {
      g.fillStyle = 'rgba(60, 177, 115, 0.2)'; // green for selected position
      g.fillRect(x + dx1 - 4, y - 5, opts.monomerTextSizeMap[colorPart].width + 8, h + 10);
      drawColor = DG.Color.toHtml(DG.Color.setAlpha(DG.Color.fromHtml(drawColor), 255));
    }
    g.fillStyle = drawColor;
    g.globalAlpha = 1.0 - opts.transparencyRate; // FIX: Use consistent transparency logic
    if (opts.drawStyle === DrawStyle.classic) {
      g.fillText(colorPart, x + dx1, y + dy);
      g.fillStyle = grayColor;
      g.fillText(grayPart, x + dx2, y + dy);
    }
    if (opts.drawStyle === DrawStyle.MSA)
      g.fillText(colorPart, x + dx1, y + dy);
    g.globalAlpha = 1.0; // Reset alpha
  }

  const placeX: number = (opts.maxWord[opts.wordIdx] ?? 0) - (opts.maxWord[0] ?? 0);
  if (opts.left || textSize > w) {
    draw(placeX, placeX + maxColorTextSize);
    return x + placeX + maxColorTextSize + grayPartSize;
  } else {
    const dx = (w - textSize) / 2;
    draw(dx, dx + maxColorTextSize);
    return x + placeX + dx + maxColorTextSize;
  }
}
