import {NUCLEOTIDES} from '../../common/model/const';
import {axolabsStyleMap} from '../../common/data-loading-utils/json-loader';

export function isOverhang(modification: string): boolean {
  const overhangSuffix = '(o)';
  return modification.endsWith(overhangSuffix);
}

export function isOneDigitNumber(n: number): boolean {
  return n >= 0 && n < 10;
}

export function getPointsToDrawStar(centerX: number, centerY: number): string {
  const outerVerticesPerStar = 5;
  const innerRadius = 3;
  const outerRadius = innerRadius * 2;
  const radiansPerVertex = Math.PI / outerVerticesPerStar;
  const radiansOffset = - radiansPerVertex / 2;
  const totalNumberOfVertices = outerVerticesPerStar * 2;

  const points = Array.from({ length: totalNumberOfVertices }, (_, i) => {
    const isOuterVertex = i % 2 === 0;
    const radius = isOuterVertex ? outerRadius : innerRadius;
    const angle = i * radiansPerVertex + radiansOffset;
    const x = centerX + Math.cos(angle) * radius;
    const y = centerY + Math.sin(angle) * radius;
    return `${x},${y}`;
  });

  return points.join(' ');
}

export function countOverhangsOnTheRightEdge(modifications: string[]): number {
  const lastIdx = modifications.length - 1;
  let count = 0;
  while (count <= lastIdx && isOverhang(modifications[count])) {
    count++;
  }
  return count === lastIdx + 1 ? 0 : count;
}

export function textWidth(text: string, fontSize: number, fontFamily: string): number {
  const canvas = document.createElement('canvas');
  const context = canvas.getContext('2d');
  if (context) {
    context.font = `${fontSize}px ${fontFamily}`;
    const metrics = context.measureText(text);
    return 2 * metrics.width;
  }
  return 0;
}

export function textInsideCircle(bases: string[], index: number): string {
  const base = bases[index];
  const isValidBase = !isOverhang(base) && NUCLEOTIDES.includes(base);
  
  return isValidBase ? base : '';
}

export function fontColorVisibleOnBackground(base: string): string {
  const RED_COEFFICIENT = 0.299;
  const GREEN_COEFFICIENT = 0.587;
  const BLUE_COEFFICIENT = 0.114;
  const LUMINANCE_THRESHOLD = 186;
  const DARK_COLOR = '#333333';
  const LIGHT_COLOR = '#ffffff';

  const styleMap = axolabsStyleMap;
  const baseColor = styleMap[base]?.color || '';
  
  const rgbValues = baseColor.match(/\d+/g)?.map(Number);
  if (!rgbValues || rgbValues.length < 3) {
    return LIGHT_COLOR;
  }

  const [r, g, b] = rgbValues;
  const luminance = r * RED_COEFFICIENT + g * GREEN_COEFFICIENT + b * BLUE_COEFFICIENT;
  return luminance > LUMINANCE_THRESHOLD ? DARK_COLOR : LIGHT_COLOR;
}

export function baseColor(base: string): string {
  const styleMap = axolabsStyleMap;
  return styleMap[base].color;
}

export const svg = {
  xmlns: 'http://www.w3.org/2000/svg',
  
  createSVGElement(type: string): Element {
    return document.createElementNS(this.xmlns, type);
  },

  render(width: number, height: number): Element {
    const e = this.createSVGElement('svg');
    e.setAttribute('id', 'mySvg');
    e.setAttribute('width', String(width));
    e.setAttribute('height', String(height));
    return e;
  },

  circle(x: number, y: number, radius: number, color: string): Element {
    const e = this.createSVGElement('circle');
    e.setAttribute('cx', String(x));
    e.setAttribute('cy', String(y));
    e.setAttribute('r', String(radius));
    e.setAttribute('fill', color);
    return e;
  },

  text(textContent: string, x: number, y: number, fontSize: number, color: string): Element {
    const e = this.createSVGElement('text');
    e.setAttribute('x', String(x));
    e.setAttribute('y', String(y));
    e.setAttribute('font-size', String(fontSize));
    e.setAttribute('font-weight', 'normal');
    e.setAttribute('font-family', 'Arial');
    e.setAttribute('fill', color);
    e.textContent = textContent;
    return e;
  },

  star(x: number, y: number, fill: string): Element {
    const e = this.createSVGElement('polygon');
    e.setAttribute('points', getPointsToDrawStar(x, y));
    e.setAttribute('fill', fill);
    return e;
  },
};
