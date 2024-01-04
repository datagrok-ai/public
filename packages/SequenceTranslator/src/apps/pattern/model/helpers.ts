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

// export function textWidth(text: string, font: number): number {
//   const context = document.createElement('canvas').getContext('2d');
//   context.font = String(font);
//   return 2 * context.measureText(text).width;
// }

export function textInsideCircle(bases: string[], index: number): string {
  return (isOverhang(bases[index]) || !NUCLEOTIDES.includes(bases[index])) ? '' : bases[index];
}

export function fontColorVisibleOnBackground(base: string): string {
  const AXOLABS_MAP = axolabsStyleMap;
  const rgbIntList = AXOLABS_MAP[base].color.match(/\d+/g)!.map((e) => Number(e));
  return (rgbIntList[0] * 0.299 + rgbIntList[1] * 0.587 + rgbIntList[2] * 0.114) > 186 ? '#33333' : '#ffffff';
}

export function baseColor(base: string): string {
  const AXOLABS_MAP = axolabsStyleMap;
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
