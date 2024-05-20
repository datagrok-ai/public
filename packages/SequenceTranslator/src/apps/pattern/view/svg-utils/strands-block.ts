import {STRAND, STRANDS} from '../../model/const';
import {PatternConfiguration} from '../../model/types';
import {SVG_CIRCLE_SIZES, SVG_ELEMENT_COLORS, SVG_TEXT_FONT_SIZES} from './const';
import {NonLegendBlockBase} from './svg-block-bases';
import {SVGElementFactory} from './svg-element-factory';
import {TextDimensionsCalculator} from './text-dimensions-calculator';
import {computeTextColorForNucleobaseLabel, getNucleobaseColorFromStyleMap, getNucleobaseLabelForCircle} from './utils';

const NUMERIC_LABEL_PADDING = 5;
const SENSE_STRAND_HEIGHT = SVG_TEXT_FONT_SIZES.NUCLEOBASE +
  NUMERIC_LABEL_PADDING + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;
const SENSE_STRAND_PADDING = 10;

export class StrandsBlock extends NonLegendBlockBase {
  private strands: NonLegendBlockBase[];
  private labels: NonLegendBlockBase[];
  constructor(
    svgElementFactory: SVGElementFactory,
    config: PatternConfiguration,
    yShift: number
  ) {
    super(svgElementFactory, config, yShift);
    const strandTypes = STRANDS.filter((strandType) => config.nucleotideSequences[strandType].length > 0);

    this.strands = strandTypes
      .map((strand) => new SingleStrandBlock(this.svgElementFactory, config, yShift, strand));

    this.labels = strandTypes.map(
      (strandType, idx) =>
        new StrandLabel(this.svgElementFactory, config, yShift, strandType, this.strands[idx] as SingleStrandBlock)
    );
  }

  get svgElements(): SVGElement[] {
    return this.strands.map((strand) => strand.svgElements).flat();
  }

  getContentWidth(): number {
    return Math.max(...this.labels.map((labelBlock) => labelBlock.getContentWidth()));
  }

  getContentHeight(): number {
    const result = this.strands
      .reduce((acc, strand) => acc + strand.getContentHeight(), 0);
    return result;
  }
}

class SingleStrandBlock extends NonLegendBlockBase {
  private _svgElements: SVGElement[];
  constructor(
    protected svgElementFactory: SVGElementFactory,
    protected config: PatternConfiguration,
    protected yShift: number,
    private strand: STRAND
  ) {
    super(svgElementFactory, config, yShift);
    this._svgElements = [
      this.createStrandCircles(),
      this.createPTOLinkageStars(),
    ].flat();
  }

  get svgElements(): SVGElement[] {
    return this._svgElements;
  }

  private getStrandCircleYShift(): number {
    return getStrandCircleYShift(this.strand, this.yShift);
  }

  private createStrandCircles(): SVGElement[] {
    const defaultShift = {
      x: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS + SENSE_STRAND_PADDING,
      y: this.getStrandCircleYShift()
    };

    const nucleotides = this.config.nucleotideSequences[this.strand];
    if (this.strand === STRAND.ANTISENSE) nucleotides.reverse();

    const elements = nucleotides
      .map((nucleotide, index) => this.createNucleotideElementGroup(nucleotide, index, defaultShift)).flat();

    return elements;
  }

  private createNucleotideElementGroup(
    nucleotide: string,
    index: number,
    defaultShift: {x: number, y: number}
  ): SVGElement[] {
    const circle = this.createNucleotideCircle(nucleotide, index, defaultShift);
    const numericLabel = this.config.nucleotidesWithNumericLabels.includes(nucleotide) ?
      this.createNucleotideNumericLabel(index, defaultShift) :
      null;

    return [circle, numericLabel].filter((element) => element !== null) as SVGElement[];
  }

  private createNucleotideCircle(
    nucleotide: string,
    index: number,
    defaultShift: {x: number, y: number}
  ): SVGElement[] {
    const color = getNucleobaseColorFromStyleMap(nucleotide);
    const centerPosition = {...defaultShift, x: defaultShift.x + index * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER};

    const circle = this.svgElementFactory
      .createCircleElement(centerPosition, SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS, color);

    const nonModifiedNucleotideLabel = this.createNucleotideLabel();

    return [circle];
  }

  private getNucleotideCircleCenterPosition(
    index: number,
    defaultShift: {x: number, y: number}
  ): {x: number, y: number} {
    return {
      x: defaultShift.x + index * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER,
      y: defaultShift.y
    };
  }

  private createNucleotideLabel(
    nucleobase: string
  ): SVGElement {
    const text = getNucleobaseLabelForCircle(nucleobase);
    const color = computeTextColorForNucleobaseLabel(nucleobase);
  }

  private getNumericLabelYShift(
    defaultShift: {x: number, y: number}
  ): number {
    return this.strand === STRAND.SENSE ?
      defaultShift.y - (SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS + NUMERIC_LABEL_PADDING) :
      defaultShift.y + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS + NUMERIC_LABEL_PADDING + SVG_TEXT_FONT_SIZES.COMMENT;
  }

  private createNucleotideNumericLabel(
    index: number,
    defaultShift: {x: number, y: number}
  ): SVGElement {
    const label = (index + 1).toString();
    const width = TextDimensionsCalculator.getTextDimensions(label, SVG_TEXT_FONT_SIZES.COMMENT).width;

    const position = {
      x: defaultShift.x + index * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER - width / 2,
      y: this.getNumericLabelYShift(defaultShift)
    };

    return this.svgElementFactory.createTextElement(
      label,
      position,
      SVG_TEXT_FONT_SIZES.COMMENT,
      SVG_ELEMENT_COLORS.TEXT
    );
  }

  private createPTOLinkageStars(): SVGElement[] {
    const ptoFlags = this.config.phosphorothioateLinkageFlags[this.strand];
    if (this.strand === STRAND.ANTISENSE) ptoFlags.reverse();

    const yShift = this.getStrandCircleYShift() + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS * 0.8;

    const elements = ptoFlags
      .map((ptoFlag, index) => {
        if (!ptoFlag) return null;

        const centerPosition = {
          x: SENSE_STRAND_PADDING + index * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER,
          y: yShift
        };

        const color = SVG_ELEMENT_COLORS.LINKAGE_STAR;
        return this.svgElementFactory.createStarElement(centerPosition, color);
      })
      .filter((element) => element !== null) as SVGElement[];

    return elements;
  }

  getContentWidth(): number {
    const numberOfMonomers = this.config.nucleotideSequences[this.strand].length;
    return numberOfMonomers * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER + 2 * SENSE_STRAND_PADDING;
  }

  getContentHeight(): number {
    return SENSE_STRAND_HEIGHT + SENSE_STRAND_PADDING;
  }
}

class StrandLabel extends NonLegendBlockBase {
  constructor(
    protected svgElementFactory: SVGElementFactory,
    protected config: PatternConfiguration,
    protected yShift: number,
    private strand: STRAND,
    private strandSvgWrapper: SingleStrandBlock
  ) {
    super(svgElementFactory, config, yShift);
    this.strandSvgWrapper.shiftElements({x: 50, y: 0});
  }

  get svgElements(): SVGElement[] {
    return [];
  }

  getContentWidth(): number {
    return this.strandSvgWrapper.getContentWidth();
  }

  getContentHeight(): number {
    return this.strandSvgWrapper.getContentHeight();
  }
}

function getStrandCircleYShift(
  strand: STRAND,
  defaultYShift: number
): number {
  return strand === STRAND.SENSE ?
    defaultYShift + NUMERIC_LABEL_PADDING + SVG_TEXT_FONT_SIZES.COMMENT + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS :
    defaultYShift + SENSE_STRAND_HEIGHT + SENSE_STRAND_PADDING + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
}
