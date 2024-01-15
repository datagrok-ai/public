import {NUCLEOTIDES} from '../../../common/model/const';
import {AXOLABS_STYLE_MAP as styleMap} from '../../../common/data-loader/json-loader';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../../model/helpers';
import { STRAND, STRANDS, TERMINUS, TERMINI } from '../../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../../model/types';
import {STRAND_END, STRAND_ENDS, LUMINANCE_COEFFICIENTS, TEXT_COLOR, SVG_CIRCLE_SIZES, SVG_TEXT_FONT_SIZES, SVG_ELEMENT_COLORS, STRAND_END_LABEL_TEXT, NUMERIC_LABEL_POSITION_OFFSET, DEFAULT_FONT_FAMILY, Y_POSITIONS_FOR_STRAND_ELEMENTS, STRAND_TO_END_TERMINUS} from './const';
import {
  computeCommentYPosition,
  computeLegendCircleYPosition,
  computeLegendTextYPosition,
  computeTotalSVGHeight,
  isSingleDigitNumber,
  countOverhangNucleotidesAtStrandEnd,
  computeTextWidthInPixels,
  getNucleobaseLabelForCircle,
  computeTextColorForNucleobaseLabel,
  getNucleobaseColorFromStyleMap,
  computeMaxWidthForStrandEnd,
} from './utils';

export class SVGRenderer {
  // private svgElementFactory = new SVGElementFactory();
  constructor(private config: PatternConfiguration) { }

  renderNucleotidePattern(): Element {
    let {
      patternName,
      isAntisenseStrandActive,
      // todo: rename
      bases: nucleobases,
      phosphorothioateLinkages,
      terminalModifications,
      comment,
      modificationsWithNumericLabels
    } = this.config;

    // todo: process copies
    nucleobases[STRAND.SENSE].reverse();
    phosphorothioateLinkages[STRAND.SENSE].reverse();

    const distinctNucleobaseTypes = [...new Set(
      nucleobases[STRAND.SENSE].concat(
        isAntisenseStrandActive ? nucleobases[STRAND.ANTISENSE] : []
      )
    )];

    const isPhosphorothioateLinkageActive = [
      ...phosphorothioateLinkages[STRAND.SENSE],
      ...(isAntisenseStrandActive ? phosphorothioateLinkages[STRAND.ANTISENSE] : [])
    ].some(linkage => linkage);

    const svgElementFactory = new SVGElementFactory();
    const svgDimensionsManager = new SVGDimensionsManager(this.config, distinctNucleobaseTypes, isPhosphorothioateLinkageActive);

    const canvasWidth = svgDimensionsManager.getCanvasWidth();
    const canvasHeight = svgDimensionsManager.getCanvasHeight();
    const image = svgElementFactory.createCanvas(canvasWidth, canvasHeight);

    const xPositionOfStrandLabels = svgDimensionsManager.getXPositionOfStrandLabels();

    const commentLabel = svgElementFactory.createTextElement(comment, xPositionOfStrandLabels[STRAND_END.LEFT], computeCommentYPosition(isAntisenseStrandActive), SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);

    const starElementLabel = isPhosphorothioateLinkageActive
      ? svgElementFactory.createStarElement(SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS, computeLegendCircleYPosition(isAntisenseStrandActive), SVG_ELEMENT_COLORS.LINKAGE_STAR)
      : null;

    const psLinkageLabel = isPhosphorothioateLinkageActive ? svgElementFactory.createTextElement('ps linkage', 2 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS - 8, computeLegendTextYPosition(isAntisenseStrandActive), SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT) : null;

    const strandEndLabels = new StrandEndLabelsManager(this.config, svgElementFactory, xPositionOfStrandLabels).getStrandLabels();
    const xPositionOfTerminusModifications = svgDimensionsManager.getXPositionOfTerminusModifications();
    const labelsTerminusModification = new TerminusModificationLabelsManager(this.config, svgElementFactory, xPositionOfTerminusModifications).getTerminusModificationLabels();

    const svgElements = [
      ...STRANDS.map(strand => Object.values(strandEndLabels[strand])),
      ...STRANDS.map(strand => Object.values(labelsTerminusModification[strand])),
    ].flat()
      .concat(commentLabel, starElementLabel, psLinkageLabel)
      .filter(element => element !== null) as SVGElement[];

    image.append(...svgElements);

    const numberOfNucleotides = STRANDS.reduce((acc, strand) => {
      acc[strand] = nucleobases[strand].filter(value => !isOverhangNucleotide(value)).length;
      return acc;
    }, {} as Record<typeof STRANDS[number], number>);

    const nucleotideCounter = {
      [STRAND.SENSE]: numberOfNucleotides[STRAND.SENSE],
      [STRAND.ANTISENSE]: numberOfNucleotides[STRAND.ANTISENSE]
    };

    function createLinkageStar(strand: StrandType, index: number): SVGElement | null {
      const isActive = phosphorothioateLinkages[strand][index];
      if (!isActive)
        return null;

      const xPosition = svgDimensionsManager.getXPositionOfLinkageStar(index, strand);
      const yPosition = svgDimensionsManager.getYPositionOfLinkageStar(strand);
      const color = SVG_ELEMENT_COLORS.LINKAGE_STAR;
      const starElement = svgElementFactory.createStarElement(xPosition, yPosition, color);

      return starElement;
    }

    function createStrandSVGElements(index: number, strand: STRAND): SVGElement[] {
      const nucleobase = nucleobases[strand][index];
      const isOverhang = isOverhangNucleotide(nucleobase);
      const yPositions = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand];

      if (!isOverhang)
        nucleotideCounter[strand]--;

      const nucleobaseIndex = strand === STRAND.SENSE ?
        index :
        nucleobases[strand].length - index;
      const displayedNumericLabelIndex = strand === STRAND.SENSE ?
        nucleotideCounter[strand] + 1 :
        numberOfNucleotides[strand] - nucleotideCounter[strand];
      const numericLabelOffset = svgDimensionsManager.computeNumericLabelXOffset(
        nucleobases[strand],
        nucleobaseIndex,
        displayedNumericLabelIndex,
      );
      const numericLabelXPosition = svgDimensionsManager.computeNucleobaseCircleXPosition(index, strand) + numericLabelOffset;

      const nucleotideNumberLabel = (!isOverhang && modificationsWithNumericLabels.includes(nucleobase))
        ? String(displayedNumericLabelIndex)
        : '';

      const nucleotideNumberTextElement = svgElementFactory.createTextElement(
        nucleotideNumberLabel, numericLabelXPosition, yPositions.NUMERIC_LABEL, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT
      );

      const nucleobaseCircleXPosition = svgDimensionsManager.computeNucleobaseCircleXPosition(index, strand);
      const nucleobaseCircleElement = svgElementFactory.createCircleElement(
        nucleobaseCircleXPosition, yPositions.NUCLEOBASE_CIRCLE, SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS, getNucleobaseColorFromStyleMap(nucleobase)
      );

      const nucleobaseLabelTextElement = svgElementFactory.createTextElement(
        getNucleobaseLabelForCircle(nucleobase), numericLabelXPosition, yPositions.NUCLEOBASE_LABEL, SVG_TEXT_FONT_SIZES.NUCLEOBASE, computeTextColorForNucleobaseLabel(nucleobase)
      );

      const phosphorothioateLinkageStarElement = createLinkageStar(strand, index);

      const lastNucleotideIndex = nucleobases[strand].length;
      const lastStar = createLinkageStar(strand, lastNucleotideIndex);

      const nucleotideSvgElements = [
        nucleotideNumberTextElement,
        nucleobaseCircleElement,
        nucleobaseLabelTextElement,
        phosphorothioateLinkageStarElement,
        lastStar,
      ].filter((element) => element !== null) as SVGElement[];

      return nucleotideSvgElements;
    }

    function appendStrandElements(strandType: STRAND) {
      nucleobases[strandType].forEach((_, index) => {
        const elements = createStrandSVGElements(index, strandType);
        image.append(...elements);
      });
    }

    appendStrandElements(STRAND.SENSE);

    if (isAntisenseStrandActive) {
      appendStrandElements(STRAND.ANTISENSE);
    }

    function createTitleElement(patternName: string,
      numberOfNucleotides: Record<typeof STRANDS[number], number>,
      isAntisenseStrandActive: boolean
    ) {
      const titleText = `${patternName} for ${numberOfNucleotides[STRAND.SENSE]}${isAntisenseStrandActive ? `/${numberOfNucleotides[STRAND.ANTISENSE]}` : ''}mer`;
      return svgElementFactory.createTextElement(
        titleText,
        SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
        SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
        SVG_TEXT_FONT_SIZES.NUCLEOBASE,
        SVG_ELEMENT_COLORS.TITLE_TEXT
      );
    }
    const titleElement = createTitleElement(patternName, numberOfNucleotides, isAntisenseStrandActive);

    image.append(titleElement);

    const legend = new LegendManager(this.config, svgDimensionsManager, svgElementFactory).getLegendItems(distinctNucleobaseTypes);
    image.append(...legend);

    return image;
  }
}

class StrandEndLabelsManager {
  constructor(
    private config: PatternConfiguration,
    private svgElementFactory: SVGElementFactory,
    private xPositionOfStrandLabels: Record<typeof STRAND_ENDS[number], number>,
  ) { }

  private createLabelForStrandEnd(strand: StrandType, end: STRAND_END): SVGElement | null {
    const isLabelActive = (strand === STRAND.SENSE) || (strand === STRAND.ANTISENSE && this.config.isAntisenseStrandActive);
    if (!isLabelActive) {
      return null;
    }

    const labelText = STRAND_END_LABEL_TEXT[end][strand];
    const xPosition = this.xPositionOfStrandLabels[end];
    const yPosition = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL;

    return this.svgElementFactory.createTextElement(labelText, xPosition, yPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.TEXT);
  }

  getStrandLabels() {
    const strandEndLabels = STRANDS.reduce((acc, strand) => {
      acc[strand] = STRAND_ENDS.reduce((endAcc, end) => {
        endAcc[end] = this.createLabelForStrandEnd(strand, end);
        return endAcc;
      }, {} as Record<typeof STRAND_ENDS[number], SVGElement | null>);
      return acc;
    }, {} as Record<typeof STRANDS[number], Record<typeof STRAND_ENDS[number], SVGElement | null>>);

    return strandEndLabels;
  }
}

class TerminusModificationLabelsManager {
  constructor(
    private config: PatternConfiguration,
    private svgElementFactory: SVGElementFactory,
    private xPositionOfTerminusModifications: Record<typeof STRAND_ENDS[number], Record<typeof STRANDS[number], number>>,
  ) { }

  private createTerminusModificationLabel(strand: StrandType, terminus: TerminalType): SVGElement | null {
    if (strand === STRAND.ANTISENSE && !this.config.isAntisenseStrandActive) {
      return null;
    }

    const end = (strand === STRAND.SENSE && terminus === TERMINUS.FIVE_PRIME) ||
      (strand === STRAND.ANTISENSE && terminus === TERMINUS.THREE_PRIME) ? STRAND_END.LEFT : STRAND_END.RIGHT;

    const labelText = this.config.terminalModifications[strand][terminus];
    const xPosition = this.xPositionOfTerminusModifications[end][strand];
    const yPosition = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL;

    return this.svgElementFactory.createTextElement(labelText, xPosition, yPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.MODIFICATION_TEXT);
  }

  getTerminusModificationLabels() {
    const labelsTerminusModification = STRANDS.reduce((acc, strand) => {
      acc[strand] = TERMINI.reduce((terminiAcc, terminus) => {
        terminiAcc[terminus] = this.createTerminusModificationLabel(strand, terminus);
        return terminiAcc;
      }, {} as Record<typeof TERMINI[number], SVGElement | null>);
      return acc;
    }, {} as Record<typeof STRANDS[number], Record<typeof TERMINI[number], SVGElement | null>>);

    return labelsTerminusModification;
  }
}

class LegendManager {
  constructor(
    private config: PatternConfiguration,
    private svgDimensionsManager: SVGDimensionsManager,
    private svgElementFactory: SVGElementFactory,
  ) { }

  private createLegendCircle(nucleobaseType: string, index: number, isAntisenseStrandActive: boolean): SVGCircleElement {
    const xPosition = this.svgDimensionsManager.computeLegendCircleXPosition(index);
    const yPosition = computeLegendCircleYPosition(isAntisenseStrandActive);
    const color = getNucleobaseColorFromStyleMap(nucleobaseType);
    return this.svgElementFactory.createCircleElement(xPosition, yPosition, SVG_CIRCLE_SIZES.LEGEND_RADIUS, color);
  }

  private createLegendText(nucleobaseType: string, index: number, isAntisenseStrandActive: boolean): SVGTextElement {
    const xPosition = this.svgDimensionsManager.computeLegendCircleXPosition(index) + SVG_CIRCLE_SIZES.LEGEND_RADIUS + 4;
    const yPosition = computeLegendTextYPosition(isAntisenseStrandActive);
    return this.svgElementFactory.createTextElement(nucleobaseType, xPosition, yPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);
  }

  getLegendItems(distinctNucleobaseTypes: string[]): SVGElement[] {
    const svgElements = [] as SVGElement[];
    distinctNucleobaseTypes.forEach((nucleobaseType, index) => {
      const legendCircle = this.createLegendCircle(nucleobaseType, index, this.config.isAntisenseStrandActive);
      const legendText = this.createLegendText(nucleobaseType, index, this.config.isAntisenseStrandActive);
      // image.append(legendCircle, legendText);
      svgElements.push(legendCircle, legendText);
    });
    return svgElements;
  }
}

class SVGDimensionsManager {
  private widthOfStrandLabel: Record<typeof STRAND_ENDS[number], number>;
  private maxWidthOfTerminusLabelsByEnd: Record<typeof STRAND_ENDS[number], number>;
  private rightOverhangs: Record<typeof STRANDS[number], number>;
  private maxStrandLength: number;
  private widthOfRightOverhangs: number;
  private xPositionOfTerminusModifications: Record<typeof STRAND_ENDS[number], Record<typeof STRANDS[number], number>>;

  constructor(
    private config: PatternConfiguration,
    // private svgElementFactory: SVGElementFactory,
    private distinctNucleobaseTypes: string[],
    private isPhosphorothioateLinkageActive: boolean,
  ) {
    this.widthOfStrandLabel = this.getWidthOfStrandLabel();
    this.maxWidthOfTerminusLabelsByEnd = this.getMaxWidthOfTerminusLabels();


    this.rightOverhangs = STRANDS.reduce((acc, strand) => {
      acc[strand] = countOverhangNucleotidesAtStrandEnd(this.config.bases[strand]);
      return acc;
    }, {} as Record<typeof STRANDS[number], number>);

    this.maxStrandLength = Math.max(...STRANDS.map(strand => this.config.bases[strand].length - this.rightOverhangs[strand]));

    this.widthOfRightOverhangs = Math.max(this.rightOverhangs[STRAND.SENSE], this.rightOverhangs[STRAND.ANTISENSE]);

    this.xPositionOfTerminusModifications = this.computeXPositionOfTerminusModifications();
  }

  // todo: more descriptive name
  getCanvasWidth(): number {
    const widthOfNucleobases = SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER * (this.maxStrandLength + this.widthOfRightOverhangs);

    const canvasWidth = STRAND_ENDS.reduce((acc, end) => {
      acc += this.widthOfStrandLabel[end] + this.maxWidthOfTerminusLabelsByEnd[end];
      return acc;
    }, 0) + widthOfNucleobases + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;

    return canvasWidth;
  }

  // todo: types to aliases
  private getWidthOfStrandLabel(): Record<typeof STRAND_ENDS[number], number> {
    const widthOfStrandLabel = Object.fromEntries(
      STRAND_ENDS.map(
        strandEnd => [strandEnd, computeMaxWidthForStrandEnd(strandEnd)]
      )
    ) as Record<typeof STRAND_ENDS[number], number>;

    return widthOfStrandLabel;
  }

  private getMaxWidthOfTerminusLabels(): Record<typeof STRAND_ENDS[number], number> {
    const maxWidthOfTerminusLabelsByEnd = STRAND_ENDS.reduce((acc, end) => {
      acc[end] = this.computeMaxWidthOfTerminalLabels(end);
      return acc;
    }, {} as Record<typeof STRAND_ENDS[number], number>);

    return maxWidthOfTerminusLabelsByEnd;
  }

  private computeMaxWidthOfTerminalLabels(end: typeof STRAND_ENDS[number]) {
    return Math.max(
      ...STRANDS.map(strand =>
        computeTextWidthInPixels(
          this.config.terminalModifications[strand][STRAND_TO_END_TERMINUS[strand][end]],
          SVG_TEXT_FONT_SIZES.NUCLEOBASE,
          DEFAULT_FONT_FAMILY
        )
      )
    );
  };


  computeLegendCircleXPosition(index: number): number {
    const legendStartIndex = this.isPhosphorothioateLinkageActive ? 1 : 0;
    const totalPositions = this.distinctNucleobaseTypes.length + legendStartIndex;
    const width = this.getCanvasWidth();
    const spacingUnit = width / totalPositions;
    const position = (index + legendStartIndex) * spacingUnit;
    const adjustedPosition = position + SVG_CIRCLE_SIZES.LEGEND_RADIUS;
    return Math.round(adjustedPosition);
  }

  // todo: cleanup legacy logic
  computeNucleobaseCircleXPosition(index: number, strand: STRAND): number {
    const rightOverhangs = this.rightOverhangs[strand];
    return this._computeNucleobaseCircleXPosition(index, rightOverhangs);
  }

  private _computeNucleobaseCircleXPosition(index: number, rightOverhangs: number): number {
    const rightModificationOffset = this.maxWidthOfTerminusLabelsByEnd[STRAND_END.RIGHT];
    const positionalIndex = this.maxStrandLength - index + rightOverhangs + 1;
    const xPosition = positionalIndex * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;
    const finalPosition = rightModificationOffset + xPosition;
    return finalPosition;
  }

  private computeXPositionOfTerminusModifications(): Record<typeof STRAND_ENDS[number], Record<typeof STRANDS[number], number>> {
    const xPositionOfTerminusModifications = Object.fromEntries(
      STRAND_ENDS.map(end => [
        end,
        Object.fromEntries(
          STRANDS.map(strand => [
            strand,
            end === STRAND_END.LEFT
            ? this.widthOfStrandLabel[STRAND_END.LEFT] - 5
            : this.rightOverhangs[strand] * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER + this._computeNucleobaseCircleXPosition(-0.5, 0)
          ])
        )
      ])
    ) as Record<typeof STRAND_ENDS[number], Record<typeof STRANDS[number], number>>;

    return xPositionOfTerminusModifications;
  }

  getXPositionOfTerminusModifications(): Record<typeof STRAND_ENDS[number], Record<typeof STRANDS[number], number>> {
    return this.xPositionOfTerminusModifications;
  }

  getXPositionOfStrandLabels(): Record<typeof STRAND_ENDS[number], number> {
    const maxRightTerminusModificationShift = Math.max(
      ...STRANDS.map(strand => this.xPositionOfTerminusModifications[STRAND_END.RIGHT][strand])
    );

    const rightEndXPosition = maxRightTerminusModificationShift
      + this.maxWidthOfTerminusLabelsByEnd[STRAND_END.LEFT]
      + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER * this.widthOfRightOverhangs;

    const xPositionOfStrandLabels = {
      [STRAND_END.LEFT]: 0,
      [STRAND_END.RIGHT]: rightEndXPosition,
    };

    return xPositionOfStrandLabels;
  }

  getCanvasHeight(): number {
    return computeTotalSVGHeight(this.config.isAntisenseStrandActive);
  }

  computeNumericLabelXOffset(nucleobases: string[], nucleobaseIndex: number, nucleotideNumericLabel: number): number {
    const criterion = isSingleDigitNumber(nucleotideNumericLabel) || NUCLEOTIDES.includes(nucleobases[nucleobaseIndex]);
    const shiftAmount = criterion
      ? NUMERIC_LABEL_POSITION_OFFSET.ONE_DIGIT
      : NUMERIC_LABEL_POSITION_OFFSET.TWO_DIGIT;
    return shiftAmount;
  }

  // todo: remove legacy division of x-y coordinates throughout the code
  getXPositionOfLinkageStar(index: number, strand: STRAND): number {
    return this._computeNucleobaseCircleXPosition(index, this.rightOverhangs[strand]) + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
  }

  getYPositionOfLinkageStar(strand: STRAND): number {
    return Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL + SVG_CIRCLE_SIZES.LINKAGE_STAR_RADIUS;
  }
}
