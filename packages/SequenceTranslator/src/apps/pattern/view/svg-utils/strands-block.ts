import * as DG from 'datagrok-api/dg';
import {NUCLEOTIDES} from '../../../common/model/const';
import {STRAND, STRANDS, TERMINUS} from '../../model/const';
import {PatternConfiguration} from '../../model/types';
import {isOverhangNucleotide} from '../../model/utils';
import {LEGEND_PADDING, SVG_CIRCLE_SIZES, SVG_ELEMENT_COLORS, SVG_TEXT_FONT_SIZES} from './const';
import {SVGBlockBase} from './svg-block-base';
import {SVGElementFactory} from './svg-element-factory';
import {TextDimensionsCalculator} from './text-dimensions-calculator';
import {computeTextColorForNucleobaseLabel, getNucleobaseColorFromStyleMap, getNucleobaseLabelForCircle} from './utils';
import { EventBus } from '../../model/event-bus';

const NUMERIC_LABEL_PADDING = 5;
const SENSE_STRAND_HEIGHT = SVG_TEXT_FONT_SIZES.NUCLEOBASE +
  NUMERIC_LABEL_PADDING + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;
const SENSE_STRAND_PADDING = 10;
const LEFT_LABEL_WIDTH = 55;
const RIGHT_LABEL_WIDTH = 20;
const MODIFICATION_LABEL_WIDTH = 50;
export const SENSE_STRAND_HORIZONTAL_SHIFT = SENSE_STRAND_PADDING + LEFT_LABEL_WIDTH;

export class StrandsBlock extends SVGBlockBase {
  private strands: SVGBlockBase[];
  private labels: SVGBlockBase[];
  constructor(
    svgElementFactory: SVGElementFactory,
    config: PatternConfiguration,
    yShift: number,
    eventBus: EventBus,
  ) {
    super(svgElementFactory, config, yShift);
    const strandTypes = STRANDS.filter((strandType) => config.nucleotideSequences[strandType].length > 0);

    yShift += MODIFICATION_LABEL_WIDTH;
    this.strands = strandTypes
      .map((strand) => new SingleStrandBlock(this.svgElementFactory, config, yShift, strand, eventBus));

    this.labels = strandTypes.map(
      (strandType, idx) =>
        new StrandLabel(this.svgElementFactory, config, yShift, strandType, this.strands[idx] as SingleStrandBlock)
    );
  }

  get svgElements(): SVGElement[] {
    const elements = [
      ...this.strands,
      ...this.labels
    ].map((block) => block.svgElements).flat();
    return elements;
  }

  getContentWidth(): number {
    return Math.max(...this.labels.map((labelBlock) => labelBlock.getContentWidth()));
  }

  getContentHeight(): number {
    const result = this.strands
      .reduce((acc, strand) => acc + strand.getContentHeight(), 0);
    return result + LEGEND_PADDING;
  }
}

class SingleStrandBlock extends SVGBlockBase {
  private _svgElements: SVGElement[];
  private nucleotideNumericLabels: (number | null)[];
  constructor(
    protected svgElementFactory: SVGElementFactory,
    protected config: PatternConfiguration,
    protected yShift: number,
    private strand: STRAND,
    private eventBus: EventBus,
  ) {
    super(svgElementFactory, config, yShift);

    // WARNING: should be computed before creating circles
    this.nucleotideNumericLabels = this.computeNucleotideNumericLabels();

    if (this.strand === STRAND.ANTISENSE) {
      this.config.phosphorothioateLinkageFlags[this.strand].reverse();
      this.config.nucleotideSequences[this.strand].reverse();
      this.nucleotideNumericLabels.reverse();
    }

    this._svgElements = [
      this.createStrandCircles(),
      this.createPTOLinkageStars(),
    ].flat();
  }

  private computeNucleotideNumericLabels(): (number | null)[] {
    let index = 0;
    const nucleotides = this.config.nucleotideSequences[this.strand];
    const indices = nucleotides.map((nucleotide) => {
      if (isOverhangNucleotide(nucleotide)) return null;
      index++;
      return index;
    });
    return indices;
  }

  get svgElements(): SVGElement[] {
    if (this.strand === STRAND.ANTISENSE && !this.config.isAntisenseStrandIncluded)
      return [];
    return this._svgElements;
  }

  private getStrandCircleYShift(): number {
    return getStrandCircleYShift(this.strand, this.yShift);
  }

  private createStrandCircles(): SVGElement[] {
    const defaultShift = {
      x: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS + SENSE_STRAND_HORIZONTAL_SHIFT,
      y: this.getStrandCircleYShift()
    };

    const nucleotides = this.config.nucleotideSequences[this.strand];

    const elements = nucleotides
      .map((nucleotide, index) => this.createNucleotideElementGroup(nucleotide, index, defaultShift)).flat();

    return elements;
  }

  private createNucleotideElementGroup(
    nucleotide: string,
    index: number,
    defaultShift: {x: number, y: number}
  ): SVGElement[] {
    const circleElements = this.createNucleotideCircleElements(nucleotide, index, defaultShift);
    const numericLabel = this.config.nucleotidesWithNumericLabels.includes(nucleotide) ?
      this.createNucleotideNumericLabel(index, defaultShift) : null;

    const modificationLabel = this.config.modificationLabelsVisible ?
      this.createNucleotideModificationLabel(nucleotide, index, defaultShift) : null;
    
    return [...circleElements, numericLabel, modificationLabel].filter((element) => element !== null) as SVGElement[];
  }

  private createNucleotideCircleElements(
    nucleotide: string,
    index: number,
    defaultShift: {x: number, y: number}
  ): (SVGElement | null)[] {
    const color = getNucleobaseColorFromStyleMap(nucleotide);
    const centerPosition = {...defaultShift, x: defaultShift.x + index * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER};

    const circle = this.svgElementFactory
      .createCircleElement(centerPosition, SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS, color, 'pointer');

    const strandIdx = this.strand === STRAND.ANTISENSE ?
      this.config.nucleotideSequences[this.strand].length - 1 - index : index;
      circle.onclick = (e) => {
      e.stopPropagation();
      e.preventDefault();
      DG.Menu.popup()
        .item('Remove', () => {
          this.eventBus.removeNucleotide(this.strand, strandIdx);
        })
        .item('Add', () => {
          this.eventBus.addNucleotide(this.strand, strandIdx);
        })
        .item('Edit', () => { })
        .show();
    };

    const nonModifiedNucleotideLetterLabel = this.createNucleotideLetterLabel(index, defaultShift, nucleotide);

    return [circle, nonModifiedNucleotideLetterLabel];
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

  private createNucleotideLetterLabel(
    index: number,
    defaultShift: {x: number, y: number},
    nucleobase: string
  ): SVGElement | null {
    if (!NUCLEOTIDES.includes(nucleobase))
      return null;

    const text = getNucleobaseLabelForCircle(nucleobase);
    const color = computeTextColorForNucleobaseLabel(nucleobase);
    // position at the very center of the circle
    const position = this.getPositionForNucleotideLabel(index, defaultShift);
    return this.svgElementFactory.createTextElement(
      text,
      position,
      SVG_TEXT_FONT_SIZES.NUCLEOBASE,
      color,
      'normal',
      'default'
    );
  }

  /** Returns the position for the letter with its center being at the center of the circle */
  private getPositionForNucleotideLabel(
    index: number,
    defaultShift: {x: number, y: number}
  ): {x: number, y: number} {
    const circleCenter = this.getNucleotideCircleCenterPosition(index, defaultShift);
    const textDimensions = TextDimensionsCalculator.getTextDimensions('A', SVG_TEXT_FONT_SIZES.NUCLEOBASE);
    return {
      x: circleCenter.x - textDimensions.width / 2,
      // the coefficient 1/3 is fine-tuned to make the text look centered
      y: circleCenter.y + textDimensions.height / 3
    };
  }

  private getNumericLabelYShift(
    defaultShift: {x: number, y: number}
  ): number {
    return this.strand === STRAND.SENSE ?
      defaultShift.y - (SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS + NUMERIC_LABEL_PADDING) :
      defaultShift.y + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS + NUMERIC_LABEL_PADDING + SVG_TEXT_FONT_SIZES.COMMENT;
  }

  private getModificationLabelYShift(
    defaultShift: {x: number, y: number},
    width: number
  ): number {
    return this.strand === STRAND.SENSE ?
      this.getNumericLabelYShift(defaultShift) - 5 - width / 2 - NUMERIC_LABEL_PADDING :
      this.getNumericLabelYShift(defaultShift) + 5 + width / 2 + NUMERIC_LABEL_PADDING;
  }

  private createNucleotideNumericLabel(
    index: number,
    defaultShift: {x: number, y: number}
  ): SVGElement | null {
    const label = this.nucleotideNumericLabels[index];
    if (label === null) return null;

    const width = TextDimensionsCalculator.getTextDimensions(label.toString(), SVG_TEXT_FONT_SIZES.COMMENT).width;

    const position = {
      x: defaultShift.x + index * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER - width / 2,
      y: this.getNumericLabelYShift(defaultShift)
    };

    return this.svgElementFactory.createTextElement(
      label.toString(),
      position,
      SVG_TEXT_FONT_SIZES.COMMENT,
      SVG_ELEMENT_COLORS.TEXT,
      'normal',
      'default'
    );
  }

  private createNucleotideModificationLabel(
    label: string,
    index: number,
    defaultShift: {x: number, y: number}
  ): SVGElement | null {
    if (label === null) return null;

    //warning! width is calculated before rotating 
    const width = TextDimensionsCalculator.getTextDimensions(label.toString(), SVG_TEXT_FONT_SIZES.MODIFICATION_LABEL).width;


    let finalLabel = label;
    let updatedWidth = width;
    const widthProportion = width / MODIFICATION_LABEL_WIDTH;
    if (widthProportion > 1) {
      const symbols = Math.ceil(label.length/widthProportion);
      finalLabel = `${label.substring(0, symbols - 3)}...`; 
      updatedWidth = TextDimensionsCalculator.getTextDimensions(finalLabel.toString(), SVG_TEXT_FONT_SIZES.MODIFICATION_LABEL).width;
    }

    const position = {
      x: defaultShift.x + index * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER - updatedWidth/2,
      y: this.getModificationLabelYShift(defaultShift, updatedWidth)
    };

    return this.svgElementFactory.createTextElement(
      finalLabel,
      position,
      SVG_TEXT_FONT_SIZES.MODIFICATION_LABEL,
      SVG_ELEMENT_COLORS.MODIFICATION_LABEL,
      'normal',
      'default',
      'rotate(-90deg)'
    );
  }

  private createPTOLinkageStars(): SVGElement[] {
    const ptoFlags = this.config.phosphorothioateLinkageFlags[this.strand];

    const yShift = this.getStrandCircleYShift() + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS * 0.8;

    const elements = ptoFlags
      .map((ptoFlag, index) => {

        const centerPosition = {
          x: SENSE_STRAND_HORIZONTAL_SHIFT + index * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER,
          y: yShift
        };

        const color = SVG_ELEMENT_COLORS.LINKAGE_STAR;

        const starElement = this.svgElementFactory.createStarElement(centerPosition, color, !ptoFlag ? '0.0' : '1.0', 'pointer');
        const strandIdx = this.strand === STRAND.ANTISENSE ?
          this.config.nucleotideSequences[this.strand].length - index : index;
        starElement.onclick = (e) => {
          e.stopPropagation();
          e.preventDefault();
          DG.Menu.popup()
            .item(!ptoFlag ? 'Add PTO' : 'Remove PTO', () => {
              this.eventBus.setPhosphorothioateLinkageFlag(this.strand, strandIdx, !ptoFlag);
            })
            .show();
        };
        
        return starElement;
      })
      .filter((element) => element !== null) as SVGElement[];

    return elements;
  }

  getContentWidth(): number {
    const numberOfMonomers = this.config.nucleotideSequences[this.strand].length;
    return numberOfMonomers * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;
  }

  getContentHeight(): number {
    if (this.strand === STRAND.ANTISENSE && !this.config.isAntisenseStrandIncluded)
      return SENSE_STRAND_PADDING;

    let height = SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER + NUMERIC_LABEL_PADDING * 2 +
      SVG_TEXT_FONT_SIZES.NUCLEOBASE;
    return height + MODIFICATION_LABEL_WIDTH;
  }
}

class StrandLabel extends SVGBlockBase {
  private _svgElements: SVGElement[];
  constructor(
    protected svgElementFactory: SVGElementFactory,
    protected config: PatternConfiguration,
    protected yShift: number,
    private strand: STRAND,
    private strandSvgWrapper: SingleStrandBlock
  ) {
    super(svgElementFactory, config, yShift);
    this._svgElements = this.createSVGElements();
    // this.strandSvgWrapper.shiftElements({x: this.getLeftLabelWidth(), y: 0});
  }

  private createSVGElements(): SVGElement[] {
    const elements = [
      this.createLeftLabel(),
      this.createRightLabel()
    ];
    return elements;
  }

  private getLeftLabelWidth(): number {
    return LEFT_LABEL_WIDTH;
  }

  private getRightLabelWidth(): number {
    return RIGHT_LABEL_WIDTH;
  }

  private createLeftLabel(): SVGTextElement {
    const terminus = this.strand === STRAND.SENSE ? TERMINUS.FIVE_PRIME : TERMINUS.THREE_PRIME;
    const text = `${this.strand}: ${terminus}  `;
    const textDimensions = TextDimensionsCalculator.getTextDimensions(text, SVG_TEXT_FONT_SIZES.NUCLEOBASE);
    const position = {
      x: SENSE_STRAND_PADDING,
      y: getStrandCircleYShift(this.strand, this.yShift) + textDimensions.height / 3
    };

    return this.svgElementFactory.createTextElement(
      text,
      position,
      SVG_TEXT_FONT_SIZES.NUCLEOBASE,
      SVG_ELEMENT_COLORS.TEXT,
      'normal',
      'default'
    );
  }

  private createRightLabel(): SVGTextElement {
    const terminus = this.strand === STRAND.SENSE ? TERMINUS.THREE_PRIME : TERMINUS.FIVE_PRIME;
    const text = `  ${terminus}`;
    const textDimensions = TextDimensionsCalculator.getTextDimensions(text, SVG_TEXT_FONT_SIZES.NUCLEOBASE);
    const position = {
      x: SENSE_STRAND_HORIZONTAL_SHIFT + this.strandSvgWrapper.getContentWidth() + 5,
      y: getStrandCircleYShift(this.strand, this.yShift) + textDimensions.height / 3
    };

    return this.svgElementFactory.createTextElement(
      text,
      position,
      SVG_TEXT_FONT_SIZES.NUCLEOBASE,
      SVG_ELEMENT_COLORS.TEXT,
      'normal',
      'default'
    );
  }

  get svgElements(): SVGElement[] {
    if (this.strand === STRAND.ANTISENSE && !this.config.isAntisenseStrandIncluded)
      return [];
    return this._svgElements;
  }

  getContentWidth(): number {
    return this.strandSvgWrapper.getContentWidth() + this.getLeftLabelWidth() + this.getRightLabelWidth() + SENSE_STRAND_PADDING;
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
