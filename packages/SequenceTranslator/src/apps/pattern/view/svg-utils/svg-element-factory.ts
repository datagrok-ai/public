import {Position} from '../types';

export class SVGElementFactory {
  private readonly xmlNamespace = 'http://www.w3.org/2000/svg';

  private createElement(elementType: string): Element {
    return document.createElementNS(this.xmlNamespace, elementType);
  }

  private setAttributes(targetElement: Element, attributes: { [key: string]: any }): void {
    Object.entries(attributes).forEach(([key, value]) => {
      targetElement.setAttribute(key, String(value));
    });
  }

  public createCanvas(width: number, height: number): SVGElement {
    const svgElement = this.createElement('svg') as SVGElement;
    this.setAttributes(svgElement, {
      id: 'mySvg',
      width,
      height,
    });
    return svgElement;
  }

  public createCircleElement(centerPosition: Position, radius: number, color: string): SVGCircleElement {
    const circle = this.createElement('circle') as SVGCircleElement;
    this.setAttributes(circle, {
      cx: centerPosition.x,
      cy: centerPosition.y,
      r: radius,
      fill: color,
    });
    return circle;
  }

  public createTextElement(textContent: string, position: Position, fontSize: number, color: string): SVGTextElement {
    const textElement = this.createElement('text') as SVGTextElement;
    this.setAttributes(textElement, {
      'x': position.x,
      'y': position.y,
      'font-size': fontSize,
      'font-weight': 'normal',
      'font-family': 'Arial',
      'fill': color,
    });
    textElement.textContent = textContent;
    return textElement;
  }

  public createStarElement(centerPosition: Position, color: string): SVGPolygonElement {
    const star = this.createElement('polygon') as SVGPolygonElement;
    const points = this.computeStarVertexCoordinates(centerPosition);
    const pointsAttribute = points.map((point) => point.join(',')).join(' ');

    this.setAttributes(star, {
      points: pointsAttribute,
      fill: color,
    });
    return star;
  }

  private computeStarVertexCoordinates(centerPosition: Position): [number, number][] {
    const outerVerticesPerStar = 5;
    const innerVertexRadius = 3;
    const outerVertexRadius = innerVertexRadius * 2;
    const radiansPerVertex = Math.PI / outerVerticesPerStar;
    const radiansOffset = - radiansPerVertex / 2;
    const totalNumberOfVertices = outerVerticesPerStar * 2;

    const points: [number, number][] = Array.from({length: totalNumberOfVertices}, (_, i) => {
      const isOuterVertex = i % 2 === 0;
      const radius = isOuterVertex ? outerVertexRadius : innerVertexRadius;
      const angle = i * radiansPerVertex + radiansOffset;
      const x = centerPosition.x + Math.cos(angle) * radius;
      const y = centerPosition.y + Math.sin(angle) * radius;
      return [x, y];
    });

    return points;
  }
}
