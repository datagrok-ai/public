import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {monomerToShort} from './macromolecule/utils';
import {ISeqSplitted} from './macromolecule/types';

const undefinedColor = 'rgb(100,100,100)';
const grayColor = '#808080';
const blackColor = 'rgb(0,0,0)';
const monomerToShortFunction: (amino: string, maxLengthOfMonomer: number) => string = monomerToShort;

export enum TAGS {
  maxMonomerLength = '.mm.cellRenderer.maxMonomerLength',
}

export enum DrawStyle {
  MSA = 'MSA',
  classic = 'classic',
}

/** A function that prints a string aligned to left or centered.
 * @param {number}x x coordinate.
 * @param {number}y y coordinate.
 * @param {number}w Width.
 * @param {number}h Height.
 * @param {CanvasRenderingContext2D}g Canvas rendering context.
 * @param {string}s String to print.
 * @param {string}color String color.
 * @param {number}pivot Pirvot.
 * @param {boolean}left Is left aligned.
 * @param {number}transparencyRate Transparency rate where 1.0 is fully transparent
 * @param {string}separator Is separator for sequence.
 * @param {boolean}last Is checker if element last or not.
 * @param {DrawStyle}drawStyle Is draw style. MSA - for multicharSeq, classic - for other seq.
 * @param {number[]}maxWord {{[pos: number]: number}}  Max word lengths per position.
 * @param {number}wordIdx Is index of word we currently draw.
 * @param {DG.GridCell}gridCell Is grid cell.
 * @param {string[]}referenceSequence Is reference sequence for diff mode.
 * @param {number}maxLengthOfMonomer Is max length of monomer.
 * @param {{[key: string]: TextMetrics}}monomerTextSizeMap Is map of monomer text sizes.
 * @return {number} x coordinate to start printing at.*/
export function printLeftOrCentered(
  x: number, y: number, w: number, h: number,
  g: CanvasRenderingContext2D, s: string, color: string = undefinedColor,
  pivot: number = 0, left: boolean = false, transparencyRate: number = 1.0,
  separator: string = '', last: boolean = false, drawStyle: DrawStyle = DrawStyle.classic,
  maxWord: number[] = [], wordIdx: number = 0, gridCell: DG.GridCell | null = null,
  referenceSequence: string[] | null = null, maxLengthOfMonomer: number | null = null,
  monomerTextSizeMap: { [key: string]: TextMetrics } = {}
): number {
  g.textAlign = 'start';
  let colorPart = s.substring(0);
  let grayPart = last ? '' : separator;
  if (drawStyle === DrawStyle.MSA) grayPart = '';
  let colorCode = true;
  let compareWithCurrent = true;
  let highlightDifference = 'difference';
  if ((gridCell != null) && (gridCell.cell.column != null)) {
    colorCode = gridCell.cell.column.temp['color-code'] ?? true;
    compareWithCurrent = gridCell.cell.column.temp['compare-with-current'] ?? true;
    highlightDifference = gridCell.cell.column.temp['highlight-difference'] ?? 'difference';
  }

  if (referenceSequence) {
    const currentMonomerCanonical = referenceSequence[wordIdx];
    if (compareWithCurrent && (referenceSequence.length > 0) && (highlightDifference === 'difference'))
      transparencyRate = (colorPart == currentMonomerCanonical) ? 0.3 : transparencyRate;
    if (compareWithCurrent && (referenceSequence.length > 0) && (highlightDifference === 'equal'))
      transparencyRate = (colorPart != currentMonomerCanonical) ? 0.3 : transparencyRate;
  }
  if (maxLengthOfMonomer != null)
    colorPart = monomerToShortFunction(colorPart, maxLengthOfMonomer);

  const fullText = colorPart + grayPart;
  monomerTextSizeMap[fullText] ??= g.measureText(fullText);
  let textSize: any = monomerTextSizeMap[fullText];
  monomerTextSizeMap[colorPart] ??= g.measureText(colorPart);
  let maxColorTextSize = monomerTextSizeMap[colorPart].width;

  monomerTextSizeMap[grayPart] ??= g.measureText(grayPart);
  const grayPartSize = monomerTextSizeMap[grayPart].width;
  const dy = h / 2 - (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2 + 1;
  textSize = textSize.width;
  if (drawStyle === DrawStyle.MSA) {
    maxColorTextSize = maxWord[wordIdx];
    textSize = maxWord[wordIdx];
  }

  /** Draw color part at {@link dx1}, and gray part at {@link dx2}. */
  function draw(dx1: number, dx2: number): void {
    const drawColor = colorCode ? color : blackColor;
    g.fillStyle = drawColor;
    g.globalAlpha = transparencyRate;
    if (drawStyle === DrawStyle.classic) {
      g.fillText(colorPart, x + dx1, y + dy);
      g.fillStyle = grayColor;
      g.fillText(grayPart, x + dx2, y + dy);
    }
    if (drawStyle === DrawStyle.MSA)
      g.fillText(colorPart, x + dx1, y + dy);
  }

  const placeX: number = (maxWord[wordIdx] ?? 0) - (maxWord[0] ?? 0);
  if (left || textSize > w) {
    draw(placeX, placeX + maxColorTextSize);
    return x + placeX + maxColorTextSize + grayPartSize;
  } else {
    const dx = (w - textSize) / 2;
    draw(dx, dx + maxColorTextSize);
    return x + placeX + dx + maxColorTextSize;
  }
}
