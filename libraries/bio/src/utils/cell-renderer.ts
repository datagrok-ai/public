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
}();

export type PrintOptions = typeof PrintOptionsDefaults;

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
  const opts = {...PrintOptionsDefaults, ...options};

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

  let colorCode = true;
  let compareWithCurrent = true;
  let highlightDifference = 'difference';
  if ((opts.gridCell != null) && (opts.gridCell.cell.column != null)) {
    colorCode = opts.gridCell.cell.column.temp['color-code'] ?? true;
    compareWithCurrent = opts.gridCell.cell.column.temp['compare-with-current'] ?? true;
    highlightDifference = opts.gridCell.cell.column.temp['highlight-difference'] ?? 'difference';
  }

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
