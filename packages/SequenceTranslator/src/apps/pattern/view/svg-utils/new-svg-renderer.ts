import _ from 'lodash';
import {PatternConfiguration} from '../../model/types';
import {MIN_WIDTH, TITLE_SHIFT} from './const';
import {StrandsBlock} from './strands-block';
import {NonLegendBlockBase, SVGBlockRendererBase} from './svg-block-bases';
import {SVGElementFactory} from './svg-element-factory';
import {TitleBlock} from './title-block';

export class NewNucleotidePatternSVGRenderer {
  private title: TitleBlock;
  private strands: StrandsBlock;
  private legend: LegendBlock;
  private svgElementFactory = new SVGElementFactory();

  constructor(patternConfig: PatternConfiguration) {
    const config = _.cloneDeep(patternConfig);

    let defaultHeight = TITLE_SHIFT;

    this.title = new TitleBlock(this.svgElementFactory, config, defaultHeight);
    defaultHeight += this.title.getContentHeight();

    this.strands = new StrandsBlock(this.svgElementFactory, config, defaultHeight);
    defaultHeight += this.strands.getContentHeight();

    this.legend = new LegendBlock(this.svgElementFactory, config, defaultHeight);
  }

  renderPattern(): SVGElement {
    const background = this.getGrayBackground();

    const width = this.getGlobalWidth();
    const height = this.getGlobalHeight();
    const canvas = this.svgElementFactory.createCanvas(width, height);

    [this.title, this.strands].forEach((block) => block.adjustContentWithinGlobalContainer(width));

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
    const blocks = [this.title, this.strands] as NonLegendBlockBase[];
    return Math.max(...blocks.map((block) => block.getContentWidth()), MIN_WIDTH);
  }

  private getGlobalHeight(): number {
    const blocks = [this.title, this.strands, this.legend];
    const height = blocks.reduce((acc, block) => acc + block.getContentHeight(), TITLE_SHIFT);
    return height;
  }
}

class LegendBlock extends SVGBlockRendererBase {
  get svgElements(): SVGElement[] {
    return [];
  }

  setWidth(width: number): void { };

  getContentHeight(): number { return 0; };
}
