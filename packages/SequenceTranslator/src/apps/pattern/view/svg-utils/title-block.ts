import {STRAND} from '../../model/const';
import {PatternConfiguration} from '../../model/types';
import {FONT_SIZE} from './const';
import {NonLegendBlockBase} from './svg-block-bases';
import {SVGElementFactory} from './svg-element-factory';
import {TextDimensionsCalculator} from './text-dimensions-calculator';

export class TitleBlock extends NonLegendBlockBase {
  private titleText: string;
  private _svgElements: SVGElement[] = [];
  constructor(
    protected svgElementFactory: SVGElementFactory,
    protected config: PatternConfiguration,
    protected yShift: number
  ) {
    super(svgElementFactory, config, yShift);
    this.titleText = this.getTitleText();
    this._svgElements = [this.getTitle()];
  }

  get svgElements(): SVGElement[] {
    return this._svgElements;
  }

  private getTitle(): SVGElement {
    return this.svgElementFactory.createTextElement(
      this.titleText,
      {x: 0, y: this.yShift + FONT_SIZE.TITLE},
      FONT_SIZE.TITLE,
      'black'
    );
  }

  getContentWidth(): number {
    return TextDimensionsCalculator.getTextDimensions(this.titleText, FONT_SIZE.TITLE).width;
  }

  getContentHeight(): number {
    return TextDimensionsCalculator.getTextDimensions(this.titleText, FONT_SIZE.TITLE).height;
  }

  private getTitleText(): string {
    const senseStrandLength = `${this.config.nucleotideSequences[STRAND.SENSE].length}`;
    const antisenseStrandLength = this.config.isAntisenseStrandIncluded ?
      `/${this.config.nucleotideSequences[STRAND.ANTISENSE].length}` : '';
    const titleText = `${this.config.patternName} for ${senseStrandLength}${antisenseStrandLength}-mer`;

    return titleText;
  }
}

