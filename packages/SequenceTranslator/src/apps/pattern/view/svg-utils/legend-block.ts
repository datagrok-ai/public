import {PatternConfiguration} from '../../model/types';
import {LEGEND_PADDING, SVG_CIRCLE_SIZES, SVG_ELEMENT_COLORS, SVG_TEXT_FONT_SIZES} from './const';
import {SVGBlockBase} from './svg-block-base';
import {SVGElementFactory} from './svg-element-factory';
import {TextDimensionsCalculator} from './text-dimensions-calculator';
import {getNucleobaseColorFromStyleMap} from './utils';

export class LegendBlock extends SVGBlockBase {
  private _svgElements: SVGElement[] = [];
  private width: number;

  constructor(
    svgElementFactory: SVGElementFactory,
    config: PatternConfiguration,
    heightShift: number
  ) {
    super(svgElementFactory, config, heightShift);
    const {elements, width} = this.createLegendItems();
    this._svgElements = elements;
    this.width = width;
  }

  private isPhosphorothioatePresent(): boolean {
    return Object.values(this.config.phosphorothioateLinkageFlags)
      .flat().some((element) => element);
  }

  private getModificationTypesList(): string[] {
    return [...new Set(
      Object.values(this.config.nucleotideSequences).flat().sort(
        (a, b) => a.toLowerCase().localeCompare(b.toLowerCase())
      )
    )];
  }

  get svgElements(): SVGElement[] {
    return this._svgElements;
  }

  private createLegendItem(
    name: string,
    xShift: number
  ): {elements: SVGElement[], width: number} {
    const circlePosition = {
      x: xShift,
      y: this.yShift - SVG_CIRCLE_SIZES.LEGEND_RADIUS
    };

    const circle = name.includes('linkage') ?
      this.svgElementFactory.createStarElement(circlePosition, SVG_ELEMENT_COLORS.LINKAGE_STAR, '1.0') :
      this.svgElementFactory
        .createCircleElement(circlePosition, SVG_CIRCLE_SIZES.LEGEND_RADIUS, getNucleobaseColorFromStyleMap(name));
    const paddedCircleWidth = 2 * SVG_CIRCLE_SIZES.LEGEND_RADIUS;
    xShift += paddedCircleWidth;

    const textPosition = {y: this.yShift, x: xShift};

    const textElement = this.svgElementFactory
      .createTextElement(name, textPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);
    const textWidth = TextDimensionsCalculator.getTextDimensions(name, SVG_TEXT_FONT_SIZES.COMMENT).width;

    return {elements: [circle, textElement], width: paddedCircleWidth + textWidth};
  }

  private createLegendItems(): { elements: SVGElement[], width: number } {
    let xShift = LEGEND_PADDING;
    const items = [] as SVGElement[][];

    const shift = (width: number) => xShift += width + 15;

    if (this.isPhosphorothioatePresent()) {
      const {elements, width} = this.createLegendItem('PTO linkage', xShift);
      shift(width);
      items.push(elements);
    }

    const modificationTypes = this.getModificationTypesList();
    modificationTypes.forEach((name) => {
      const {elements, width} = this.createLegendItem(name, xShift);
      shift(width);
      items.push(elements);
    });

    return {elements: items.flat(), width: xShift};
  }

  getContentHeight(): number { return 20; };

  getContentWidth(): number {
    return this.width;
  }
}
