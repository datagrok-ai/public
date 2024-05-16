import {PatternConfiguration} from '../../model/types';
import {SVGElementFactory} from './svg-element-factory';

/** Horizontal block within SVG: title, strands, legend  */
export abstract class SVGBlockRendererBase {
  // protected svgElements: SVGElement[] = [];
  constructor(
    protected svgElementFactory: SVGElementFactory,
    protected config: PatternConfiguration,
    protected yShift: number
  ) {}

  abstract get svgElements(): SVGElement[];

  abstract getContentHeight(): number;

  shiftElements(shift: {x: number, y: number}): void {
    this.svgElements.forEach((element) => {
      const transform = element.getAttribute('transform') || '';
      const match = transform.match(/translate\(([^,]+),([^,]+)\)/);
      const x = match ? parseFloat(match[1]) : 0;
      const y = match ? parseFloat(match[2]) : 0;
      const newTransform = `translate(${x + shift.x},${y + shift.y})`;
      element.setAttribute('transform', `${transform} ${newTransform}`);
    });
  }
}

/**  For width defining blocks, like title and strands */
export abstract class NonLegendBlockBase extends SVGBlockRendererBase {
  adjustContentWithinGlobalContainer(globalWidth: number): void {
    const contentWidth = this.getContentWidth();
    if (contentWidth < globalWidth) {
      const shift = (globalWidth - contentWidth) / 2;
      this.shiftElements({x: shift, y: 0});
    }
  }

  abstract getContentWidth(): number;
}

