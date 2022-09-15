const undefinedColor = 'rgb(100,100,100)';
const grayColor = '#808080';

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
 * @param maxWordIdx Is index of word we currently draw.
 * @param gridCell Is grid cell, new for updating data in maxWord while rendering.
 * @return {number} x coordinate to start printing at.
 */
export function printLeftOrCentered(
  x: number, y: number, w: number, h: number,
  g: CanvasRenderingContext2D, s: string, color = undefinedColor,
  pivot: number = 0, left = false, transparencyRate: number = 1.0,
  separator: string = '', last: boolean = false, drawStyle: DrawStyle = DrawStyle.classic, maxWord: { [index: string]: number } = {}, maxWordIdx: number = 0, gridCell: any = {}): number {
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
    maxColorTextSize = maxWord[maxWordIdx];
    textSize = maxWord[maxWordIdx];
    if (maxWordIdx > (maxWord['bio-maxIndex'] ?? 0)) {
      maxWord['bio-maxIndex'] = maxWordIdx;
      gridCell.cell.column.temp = maxWord;
    }
  }

  function draw(dx1: number, dx2: number): void {
    g.fillStyle = color;
    g.globalAlpha = transparencyRate;
    if (drawStyle === DrawStyle.classic) {
      g.fillText(colorPart, x + dx1, y + dy);
      g.fillStyle = grayColor;
      g.fillText(grayPart, x + dx2, y + dy);
    }
    if (drawStyle === DrawStyle.MSA) {
      g.fillStyle = color;
      g.fillText(colorPart, x + dx1 + ((maxWord[maxWordIdx] - colorTextSize) / 2), y + dy);
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

