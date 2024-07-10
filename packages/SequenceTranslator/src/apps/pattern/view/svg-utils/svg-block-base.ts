import {PatternConfiguration} from '../../model/types';
import {SVGElementFactory} from './svg-element-factory';

/** Horizontal block within SVG: title, strands, legend  */
export abstract class SVGBlockBase {
  // protected svgElements: SVGElement[] = [];
  constructor(
    protected svgElementFactory: SVGElementFactory,
    protected config: PatternConfiguration,
    protected yShift: number
  ) {}

  abstract get svgElements(): SVGElement[];

  abstract getContentHeight(): number;

  shiftElements(shift: { x: number, y: number }): void {
    this.svgElements.forEach((element) => {
      if (element.getAttribute('x'))
        element.setAttribute('x', `${parseFloat(element.getAttribute('x')!) + shift.x}`);
      else if (element.getAttribute('cx'))
        element.setAttribute('cx', `${parseFloat(element.getAttribute('cx')!) + shift.x}`);
      else {
        const points = element.getAttribute('points');
        if (points) {
          const starCoords = points.split(',').map((it) => parseFloat(it));
          element.setAttribute('points', `${starCoords[0] + shift.x}, ${starCoords[1]}`);
        }
      }
    });
  }

  abstract getContentWidth(): number;
}
