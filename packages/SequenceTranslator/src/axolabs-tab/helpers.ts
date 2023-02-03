import {AXOLABS_MAP} from '../hardcode-to-be-eliminated/constants';
import {NUCLEOTIDES} from '../hardcode-to-be-eliminated/map';

export function isOverhang(modification: string): boolean {
  return modification.slice(-3) == '(o)';
}

export function isOneDigitNumber(n: number): boolean {
  return n < 10;
}

// https://uxdesign.cc/star-rating-make-svg-great-again-d4ce4731347e
export function getPointsToDrawStar(centerX: number, centerY: number): string {
  const innerCirclePoints = 5; // a 5 point star
  const innerRadius = 15 / innerCirclePoints;
  const innerOuterRadiusRatio = 2; // outter circle is x2 the inner
  const outerRadius = innerRadius * innerOuterRadiusRatio;
  const angle = Math.PI / innerCirclePoints;
  const angleOffsetToCenterStar = 60;
  const totalNumberOfPoints = innerCirclePoints * 2; // 10 in a 5-points star

  let points = '';
  for (let i = 0; i < totalNumberOfPoints; i++) {
    const r = (i % 2 == 0) ? outerRadius : innerRadius;
    const currentX = centerX + Math.cos(i * angle + angleOffsetToCenterStar) * r;
    const currentY = centerY + Math.sin(i * angle + angleOffsetToCenterStar) * r;
    points += `${currentX},${currentY} `;
  }
  return points;
}

export function countOverhangsOnTheRightEdge(modifications: string[]): number {
  let i = 0;
  while (i < modifications.length && isOverhang(modifications[i]))
    i++;
  return (i == modifications.length - 1) ? 0 : i;
}

export function textWidth(text: string, font: number): number {
  const context = document.createElement('canvas').getContext('2d');
  // @ts-ignore
  context.font = String(font);
  // @ts-ignore
  return 2 * context.measureText(text).width;
}

export function textInsideCircle(bases: string[], index: number): string {
  return (isOverhang(bases[index]) || !NUCLEOTIDES.includes(bases[index])) ? '' : bases[index];
}

export function fontColorVisibleOnBackground(base: string): string {
  const rgbIntList = AXOLABS_MAP[base].color.match(/\d+/g)!.map((e) => Number(e));
  return (rgbIntList[0] * 0.299 + rgbIntList[1] * 0.587 + rgbIntList[2] * 0.114) > 186 ? '#33333' : '#ffffff';
}

export function baseColor(base: string): string {
  return AXOLABS_MAP[base].color;
}

export const svg = {
  xmlns: 'http://www.w3.org/2000/svg',
  render: function(width: number, height: number): Element {
    const e = document.createElementNS(this.xmlns, 'svg');
    e.setAttribute('id', 'mySvg');
    e.setAttribute('width', String(width));
    e.setAttribute('height', String(height));
    return e;
  },
  circle: function(x: number, y: number, radius: number, color: string): Element {
    const e = document.createElementNS(this.xmlns, 'circle');
    e.setAttribute('cx', String(x));
    e.setAttribute('cy', String(y));
    e.setAttribute('r', String(radius));
    e.setAttribute('fill', color);
    return e;
  },
  text: function(text: string, x: number, y: number, fontSize: number, color: string): Element {
    const e = document.createElementNS(this.xmlns, 'text');
    e.setAttribute('x', String(x));
    e.setAttribute('y', String(y));
    e.setAttribute('font-size', String(fontSize));
    e.setAttribute('font-weight', 'normal');
    e.setAttribute('font-family', 'Arial');
    e.setAttribute('fill', color);
    e.innerHTML = text;
    return e;
  },
  star: function(x: number, y: number, fill: string): Element {
    const e = document.createElementNS(this.xmlns, 'polygon');
    e.setAttribute('points', getPointsToDrawStar(x, y));
    e.setAttribute('fill', fill);
    return e;
  },
};
