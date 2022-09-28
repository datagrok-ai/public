import * as DG from 'datagrok-api/dg';

const undefinedColor = 'rgb(100,100,100)';
const grayColor = '#808080';
const blackColor = 'rgb(0,0,0)';

export enum DrawStyle {
  MSA = 'MSA',
  classic = 'classic',
}

/**
 * A function that prints a string aligned to left or centered.
 *
 * @param {number} x x coordinate.
 * @param {number} y y coordinate.
 * @param {number} w Width.
 * @param {number} h Height.
 * @param {CanvasRenderingContext2D} g Canvas rendering context.
 * @param {string} s String to print.
 * @param {string} [color=undefinedColor] String color.
 * @param {number} [pivot=0] Pirvot.
 * @param {boolean} [left=false] Is left aligned.
 * @param {number} [transparencyRate=0.0] Transparency rate where 1.0 is fully transparent
 * @param {string} [separator=''] Is separator for sequence.
 * @param {boolean} [last=false] Is checker if element last or not.
 * @param drawStyle Is draw style. MSA - for multicharSeq, classic - for other seq.
 * @param maxWord Is array of max words for each line.
 * @param wordIdx Is index of word we currently draw.
 * @param gridCell Is grid cell, new for updating data in maxWord while rendering.
 * @param { [index: string]: number }  additionalData Is additional data for rendering.
 * @return {number} x coordinate to start printing at.
 */
export function printLeftOrCentered(
  x: number, y: number, w: number, h: number,
  g: CanvasRenderingContext2D, s: string, color = undefinedColor,
  pivot: number = 0, left = false, transparencyRate: number = 1.0,
  separator: string = '', last: boolean = false, drawStyle: DrawStyle = DrawStyle.classic, maxWord: { [index: string]: number } = {}, wordIdx: number = 0, gridCell: DG.GridCell | null = null, additionalData: { [index: string]: boolean | number | string } = {}): number {
  g.textAlign = 'start';
  const colorPart = s.substring(0);
  let grayPart = last ? '' : separator;
  if (drawStyle === DrawStyle.MSA) {
    grayPart = '';
  }

  let textSize: any = g.measureText(colorPart + grayPart);
  const indent = 5;

  let maxColorTextSize = g.measureText(colorPart).width;
  let colorTextSize = g.measureText(colorPart).width;
  const dy = (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent) / 2;
  textSize = textSize.width;
  if (drawStyle === DrawStyle.MSA) {
    maxColorTextSize = maxWord[wordIdx];
    textSize = maxWord[wordIdx];
  }
  const currentMonomer = String(additionalData['referenceSequence'])[wordIdx];
  if ((additionalData['compareWithCurrent'] == true) && (String(additionalData['referenceSequence']).length > 0)) {
      transparencyRate = (colorPart == currentMonomer) ? transparencyRate : 0.4;
  }


  function draw(dx1: number, dx2: number): void {
    const drawColor = (additionalData['colorCode'] == true) ? color : blackColor;
    g.fillStyle = drawColor;
    g.globalAlpha = transparencyRate;
    if (drawStyle === DrawStyle.classic) {
      g.fillText(colorPart, x + dx1, y + dy);
      g.fillStyle = grayColor;
      g.fillText(grayPart, x + dx2, y + dy);
    }
    if (drawStyle === DrawStyle.MSA) {
      g.fillStyle = drawColor;
      g.fillText(colorPart, x + dx1 + ((maxWord[wordIdx] - colorTextSize) / 2), y + dy);
    }
  }

  if (left || textSize > w) {
    draw(indent, indent + maxColorTextSize);
    return x + maxColorTextSize + g.measureText(grayPart).width;

  } else {
    const dx = (w - textSize) / 2;
    draw(dx, dx + maxColorTextSize);
    return x + dx + maxColorTextSize;
  }
}

