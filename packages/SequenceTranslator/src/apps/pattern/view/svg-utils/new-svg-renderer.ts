import _ from 'lodash';
import {PatternConfiguration} from '../../model/types';
import {SVG_CIRCLE_SIZES, SVG_ELEMENT_COLORS, SVG_TEXT_FONT_SIZES, TITLE_SHIFT} from './const';
import {StrandsBlock} from './strands-block';
import {SVGBlockBase} from './svg-block-bases';
import {SVGElementFactory} from './svg-element-factory';
import {TitleBlock} from './title-block';
import {TextDimensionsCalculator} from './text-dimensions-calculator';
import {getNucleobaseColorFromStyleMap} from './utils';

const LEGEND_PADDING = 10;

export class NewNucleotidePatternSVGRenderer {
  private title: TitleBlock;
  private strands: StrandsBlock;
  private legend: LegendBlock;
  private svgElementFactory = new SVGElementFactory();

  constructor(patternConfig: PatternConfiguration) {
    const config = _.cloneDeep(patternConfig);

    let heightShift = TITLE_SHIFT;

    this.title = new TitleBlock(this.svgElementFactory, config, heightShift);
    heightShift += this.title.getContentHeight();

    this.strands = new StrandsBlock(this.svgElementFactory, config, heightShift);
    heightShift += this.strands.getContentHeight() + LEGEND_PADDING;

    this.legend = new LegendBlock(this.svgElementFactory, config, heightShift);
  }

  renderPattern(): SVGElement {
    const background = this.getGrayBackground();

    const width = this.getGlobalWidth();
    const height = this.getGlobalHeight();
    const canvas = this.svgElementFactory.createCanvas(width, height);

    const elements = [this.title, this.strands, this.legend].map((block) => block.svgElements).flat();
    canvas.append(background, ...elements);

    return canvas;
  }

  private getGrayBackground(): SVGElement {
    const width = this.getGlobalWidth();
    const height = this.getGlobalHeight();
    const background = this.svgElementFactory.createRectangleElement({x: 0, y: 0}, width, height, '#C0C0C0');
    return background;
  }

  private getGlobalWidth(): number {
    const blocks = [this.title, this.strands, this.legend] as SVGBlockBase[];
    return Math.max(...blocks.map((block) => block.getContentWidth()));
  }

  private getGlobalHeight(): number {
    const blocks = [this.title, this.strands, this.legend];
    const height = blocks.reduce((acc, block) => acc + block.getContentHeight(), TITLE_SHIFT);
    return height;
  }
}

class LegendBlock extends SVGBlockBase {
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
      this.svgElementFactory.createStarElement(circlePosition, SVG_ELEMENT_COLORS.LINKAGE_STAR) :
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
