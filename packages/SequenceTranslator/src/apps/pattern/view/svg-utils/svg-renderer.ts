import {NUCLEOTIDES} from '../../../common/model/const';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../../model/helpers';
import { STRAND, STRANDS, STRAND_END, STRAND_TO_END_TERMINUS_MAP, STRAND_ENDS, TERMINUS, TERMINI } from '../../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../../model/types';
import { LUMINANCE_COEFFICIENTS, TEXT_COLOR, SVG_CIRCLE_SIZES, SVG_TEXT_FONT_SIZES, SVG_ELEMENT_COLORS, STRAND_END_LABEL_TEXT, NUMERIC_LABEL_POSITION_OFFSET, DEFAULT_FONT_FAMILY, Y_POSITIONS_FOR_STRAND_ELEMENTS} from './const';
import {Position, StrandToNumberMap, StrandEndToNumberMap, StrandEndToSVGElementsMap, TerminusToSVGElementMap} from '../types';
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

export class NucleotidePatternSVGRenderer {
  private svgElementFactory = new SVGElementFactory();
  private config: PatternConfiguration;
  private isPhosphorothioateLinkageActive: boolean;
  private distinctNucleobaseTypes: string[];
  private svgDimensionsManager: PatternSVGDimensionsCalculator;

  constructor(patternConfig: PatternConfiguration) {
    this.initializeConfig(patternConfig);

    this.distinctNucleobaseTypes = this.extractUniqueNucleobases();
    this.isPhosphorothioateLinkageActive = this.checkAnyPhosphorothioateLinkages();

    this.svgDimensionsManager = new PatternSVGDimensionsCalculator(this.config, this.distinctNucleobaseTypes, this.isPhosphorothioateLinkageActive);
  }

  renderNucleotidePattern(): Element {
    let {
      patternName,
      isAntisenseStrandIncluded,
      // todo: rename to nucleobases
      nucleotideSequences,
      phosphorothioateLinkageFlags,
      strandTerminusModifications,
      patternComment,
      nucleotidesWithNumericLabels
    } = this.config;

    const patternSVGCanvas = this.createCanvas();

    const commentLabel = this.createCommentLabel();

    const starElementLabel = this.createStarElementLabel(this.isPhosphorothioateLinkageActive);

    const phosphorothioateLinkageLabel = this.getPhosphorothioateLinkageLabel(this.isPhosphorothioateLinkageActive);

    const strandEndLabels = this.getStrandEndLabels();

    const labelsTerminusModification = this.getTerminusModificationLabels();

    const labelElements = [
      ...STRANDS.flatMap(strand => [
        ...Object.values(strandEndLabels[strand]),
        ...Object.values(labelsTerminusModification[strand])
      ]),
      commentLabel,
      starElementLabel,
      phosphorothioateLinkageLabel
    ].filter(element => element !== null) as SVGElement[];

    patternSVGCanvas.append(...labelElements);

    const numberOfNucleotides = STRANDS.reduce((acc, strand) => {
      acc[strand] = nucleotideSequences[strand].filter(value => !isOverhangNucleotide(value)).length;
      return acc;
    }, {} as StrandToNumberMap);

    const nucleotideCounter = {
      [STRAND.SENSE]: numberOfNucleotides[STRAND.SENSE],
      [STRAND.ANTISENSE]: numberOfNucleotides[STRAND.ANTISENSE]
    };

    const svgDimensionsManager = this.svgDimensionsManager;
    const svgElementFactory = this.svgElementFactory;

    function createPhosphorothioateLinkageStar(strand: StrandType, index: number): SVGElement | null {
      const isActive = phosphorothioateLinkageFlags[strand][index];
      if (!isActive)
        return null;

      const centerPosition = svgDimensionsManager.getCenterPositionOfLinkageStar(index, strand);
      const color = SVG_ELEMENT_COLORS.LINKAGE_STAR;
      const starElement = svgElementFactory.createStarElement(centerPosition, color);

      return starElement;
    }

    function generateStrandNucleotideElements(index: number, strand: STRAND): SVGElement[] {
      const nucleobase = nucleotideSequences[strand][index];
      const isOverhang = isOverhangNucleotide(nucleobase);
      const yPositions = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand];

      if (!isOverhang)
        nucleotideCounter[strand]--;

      const nucleobaseIndex = strand === STRAND.SENSE ?
        index :
        nucleotideSequences[strand].length - index;
      const displayedNumericLabelIndex = strand === STRAND.SENSE ?
        nucleotideCounter[strand] + 1 :
        numberOfNucleotides[strand] - nucleotideCounter[strand];
      const numericLabelOffset = svgDimensionsManager.computeNumericLabelXOffset(
        nucleotideSequences[strand],
        nucleobaseIndex,
        displayedNumericLabelIndex,
      );
      const numericLabelXPosition = svgDimensionsManager.computeNucleobaseCircleXPosition(index, strand) + numericLabelOffset;

      const nucleotideNumberLabel = (!isOverhang && nucleotidesWithNumericLabels.includes(nucleobase))
        ? String(displayedNumericLabelIndex)
        : '';

      const nucleotideNumberTextElementPosition = {
        x: numericLabelXPosition,
        y: yPositions.NUMERIC_LABEL,
      };
      const nucleotideNumberTextElement = svgElementFactory.createTextElement(
        nucleotideNumberLabel, nucleotideNumberTextElementPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT
      );

      const nucleobaseCirclePosition = {
        x: svgDimensionsManager.computeNucleobaseCircleXPosition(index, strand),
        y: yPositions.NUCLEOBASE_CIRCLE,
      };
      const nucleobaseCircleElement = svgElementFactory.createCircleElement(
        nucleobaseCirclePosition, SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS, getNucleobaseColorFromStyleMap(nucleobase)
      );

      const nucleobaseLabelTextPosition = {
        x: numericLabelXPosition,
        y: yPositions.NUCLEOBASE_LABEL,
      };
      const nucleobaseLabelTextElement = svgElementFactory.createTextElement(
        getNucleobaseLabelForCircle(nucleobase), nucleobaseLabelTextPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, computeTextColorForNucleobaseLabel(nucleobase)
      );

      const phosphorothioateLinkageStarElement = createPhosphorothioateLinkageStar(strand, index);

      const lastNucleotideIndex = nucleotideSequences[strand].length;
      const lastStar = createPhosphorothioateLinkageStar(strand, lastNucleotideIndex);

      const nucleotideSvgElements = [
        nucleotideNumberTextElement,
        nucleobaseCircleElement,
        nucleobaseLabelTextElement,
        phosphorothioateLinkageStarElement,
        lastStar,
      ].filter((element) => element !== null) as SVGElement[];

      return nucleotideSvgElements;
    }

    function appendStrandElementsToCanvas(strandType: STRAND) {
      nucleotideSequences[strandType].forEach((_, index) => {
        const elements = generateStrandNucleotideElements(index, strandType);
        patternSVGCanvas.append(...elements);
      });
    }

    appendStrandElementsToCanvas(STRAND.SENSE);

    if (isAntisenseStrandIncluded) {
      appendStrandElementsToCanvas(STRAND.ANTISENSE);
    }

    function createTitleElement(patternName: string,
      numberOfNucleotides: StrandToNumberMap,
      isAntisenseStrandActive: boolean
    ) {
      const titleText = `${patternName} for ${numberOfNucleotides[STRAND.SENSE]}${isAntisenseStrandActive ? `/${numberOfNucleotides[STRAND.ANTISENSE]}` : ''}mer`;
      const titleTextPosition = {
        x: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
        y: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
      };
      return svgElementFactory.createTextElement(
        titleText,
        titleTextPosition,
        SVG_TEXT_FONT_SIZES.NUCLEOBASE,
        SVG_ELEMENT_COLORS.TITLE_TEXT
      );
    }
    const titleElement = createTitleElement(patternName, numberOfNucleotides, isAntisenseStrandIncluded);

    patternSVGCanvas.append(titleElement);

    const legend = new LegendManager(this.config, svgDimensionsManager, svgElementFactory).getLegendItems(this.distinctNucleobaseTypes);
    patternSVGCanvas.append(...legend);

    return patternSVGCanvas;
  }

  private extractUniqueNucleobases(): string[] {
    const distinctNucleobaseTypes = [...new Set(
      this.config.nucleotideSequences[STRAND.SENSE].concat(
        this.config.isAntisenseStrandIncluded ? this.config.nucleotideSequences[STRAND.ANTISENSE] : []
      )
    )];

    return distinctNucleobaseTypes;
  }

  private checkAnyPhosphorothioateLinkages(): boolean {
    return [...this.config.phosphorothioateLinkageFlags[STRAND.SENSE],
      ...(this.config.isAntisenseStrandIncluded ? this.config.phosphorothioateLinkageFlags[STRAND.ANTISENSE] : [])
  ].some(linkage => linkage);
  }

  private initializeConfig(patternConfig: PatternConfiguration): void {
    // WARNING: to ensure immutability, we need to deep copy the config object
    this.config = JSON.parse(JSON.stringify(patternConfig)) as PatternConfiguration;

    this.config.nucleotideSequences[STRAND.SENSE].reverse();
    this.config.phosphorothioateLinkageFlags[STRAND.SENSE].reverse();
  }

  private createCanvas(): SVGElement {
    const canvasWidth = this.svgDimensionsManager.getCanvasWidth();
    const canvasHeight = this.svgDimensionsManager.getCanvasHeight();
    const svgCanvas = this.svgElementFactory.createCanvas(canvasWidth, canvasHeight);

    return svgCanvas;
  }

  private createCommentLabel(): SVGElement {
    const commentLabelPosition = this.svgDimensionsManager.getCommentLabelPosition();

    const commentLabel = this.svgElementFactory.createTextElement(
      this.config.patternComment,
      commentLabelPosition,
      SVG_TEXT_FONT_SIZES.COMMENT,
      SVG_ELEMENT_COLORS.TEXT
    );

    return commentLabel;
  }

  createStarElementLabel(isPhosphorothioateLinkageActive: boolean): SVGElement | null {
    const starElementLabelPosition = this.svgDimensionsManager.getStarElementLabelPosition();
    const starElementLabel = isPhosphorothioateLinkageActive
      ? this.svgElementFactory.createStarElement(starElementLabelPosition, SVG_ELEMENT_COLORS.LINKAGE_STAR)
      : null;

    return starElementLabel;
  }

  getPhosphorothioateLinkageLabel(isPhosphorothioateLinkageActive: boolean): SVGElement | null{
    const position = this.svgDimensionsManager.getPhosphorothioateLinkageLabelPosition();
    const phosphorothioateLinkageLabel = isPhosphorothioateLinkageActive
      ? this.svgElementFactory.createTextElement(
        'ps linkage',
        position,
        SVG_TEXT_FONT_SIZES.COMMENT,
        SVG_ELEMENT_COLORS.TEXT)
      : null;

    return phosphorothioateLinkageLabel;
  };

  private createLabelForStrandEnd(strand: StrandType, end: STRAND_END): SVGElement | null {
    const isLabelActive = (strand === STRAND.SENSE) || (strand === STRAND.ANTISENSE && this.config.isAntisenseStrandIncluded);
    if (!isLabelActive) {
      return null;
    }

    const labelText = STRAND_END_LABEL_TEXT[end][strand];
    const labelPosition = this.svgDimensionsManager.getStrandEndLabelPosition(strand, end);

    return this.svgElementFactory.createTextElement(labelText, labelPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.TEXT);
  }

  private getStrandEndLabels() {
    const strandEndLabels = STRANDS.reduce((acc, strand) => {
      acc[strand] = STRAND_ENDS.reduce((endAcc, end) => {
        endAcc[end] = this.createLabelForStrandEnd(strand, end);
        return endAcc;
      }, {} as StrandEndToSVGElementsMap);
      return acc;
    }, {} as Record<typeof STRANDS[number], StrandEndToSVGElementsMap>);

    return strandEndLabels;
  }

  private createTerminusModificationLabel(strand: StrandType, terminus: TerminalType): SVGElement | null {
    if (strand === STRAND.ANTISENSE && !this.config.isAntisenseStrandIncluded) {
      return null;
    }

    const end = (strand === STRAND.SENSE && terminus === TERMINUS.FIVE_PRIME) ||
      (strand === STRAND.ANTISENSE && terminus === TERMINUS.THREE_PRIME) ? STRAND_END.LEFT : STRAND_END.RIGHT;

    const labelText = this.config.strandTerminusModifications[strand][terminus];
    const labelPosition = this.svgDimensionsManager.getTerminalModificationLabelPosition(strand, end);
    return this.svgElementFactory.createTextElement(labelText, labelPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.MODIFICATION_TEXT);
  }

  private getTerminusModificationLabels() {
    const labelsTerminusModification = STRANDS.reduce((acc, strand) => {
      acc[strand] = TERMINI.reduce((terminiAcc, terminus) => {
        terminiAcc[terminus] = this.createTerminusModificationLabel(strand, terminus);
        return terminiAcc;
      }, {} as TerminusToSVGElementMap);
      return acc;
    }, {} as Record<typeof STRANDS[number], TerminusToSVGElementMap>);

    return labelsTerminusModification;
  }
}

class TerminusModificationLabelsManager {
  constructor(
    private config: PatternConfiguration,
    private svgElementFactory: SVGElementFactory,
    private xPositionOfTerminusModifications: Record<typeof STRAND_ENDS[number], StrandToNumberMap>,
  ) { }

}

class LegendManager {
  constructor(
    private config: PatternConfiguration,
    private svgDimensionsManager: PatternSVGDimensionsCalculator,
    private svgElementFactory: SVGElementFactory,
  ) { }

  private createLegendCircle(nucleobaseType: string, index: number, isAntisenseStrandActive: boolean): SVGCircleElement {
    const centerPosition = {
      x: this.svgDimensionsManager.computeLegendCircleXPosition(index),
      y: computeLegendCircleYPosition(isAntisenseStrandActive),
    };
    const radius = SVG_CIRCLE_SIZES.LEGEND_RADIUS;
    const color = getNucleobaseColorFromStyleMap(nucleobaseType);

    return this.svgElementFactory.createCircleElement(centerPosition, radius, color);
  }

  private createLegendText(nucleobaseType: string, index: number, isAntisenseStrandActive: boolean): SVGTextElement {
    const legendPosition = {
      x: this.svgDimensionsManager.computeLegendCircleXPosition(index) + SVG_CIRCLE_SIZES.LEGEND_RADIUS + 4,
      y: computeLegendTextYPosition(isAntisenseStrandActive),
    };
    return this.svgElementFactory.createTextElement(nucleobaseType, legendPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);
  }

  getLegendItems(distinctNucleobaseTypes: string[]): SVGElement[] {
    const svgElements = [] as SVGElement[];
    distinctNucleobaseTypes.forEach((nucleobaseType, index) => {
      const legendCircle = this.createLegendCircle(nucleobaseType, index, this.config.isAntisenseStrandIncluded);
      const legendText = this.createLegendText(nucleobaseType, index, this.config.isAntisenseStrandIncluded);
      svgElements.push(legendCircle, legendText);
    });
    return svgElements;
  }
}

class PatternSVGDimensionsCalculator {
  private widthOfStrandLabel: StrandEndToNumberMap;
  private maxWidthOfTerminusLabelsByEnd: StrandEndToNumberMap;
  private rightOverhangs: StrandToNumberMap;
  private maxStrandLength: number;
  private widthOfRightOverhangs: number;
  private xPositionOfTerminusModifications: Record<typeof STRAND_ENDS[number], StrandToNumberMap>;

  constructor(
    private config: PatternConfiguration,
    private distinctNucleobaseTypes: string[],
    private isPhosphorothioateLinkageActive: boolean,
  ) {
    this.widthOfStrandLabel = this.getWidthOfStrandLabel();
    this.maxWidthOfTerminusLabelsByEnd = this.getMaxWidthOfTerminusLabels();


    this.rightOverhangs = STRANDS.reduce((acc, strand) => {
      acc[strand] = countOverhangNucleotidesAtStrandEnd(this.config.nucleotideSequences[strand]);
      return acc;
    }, {} as StrandToNumberMap);

    this.maxStrandLength = Math.max(...STRANDS.map(strand => this.config.nucleotideSequences[strand].length - this.rightOverhangs[strand]));

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
  private getWidthOfStrandLabel(): StrandEndToNumberMap {
    const widthOfStrandLabel = Object.fromEntries(
      STRAND_ENDS.map(
        strandEnd => [strandEnd, computeMaxWidthForStrandEnd(strandEnd)]
      )
    ) as StrandEndToNumberMap;

    return widthOfStrandLabel;
  }

  private getMaxWidthOfTerminusLabels(): StrandEndToNumberMap {
    const maxWidthOfTerminusLabelsByEnd = STRAND_ENDS.reduce((acc, end) => {
      acc[end] = this.computeMaxWidthOfTerminalLabels(end);
      return acc;
    }, {} as StrandEndToNumberMap);

    return maxWidthOfTerminusLabelsByEnd;
  }

  private computeMaxWidthOfTerminalLabels(end: typeof STRAND_ENDS[number]) {
    return Math.max(
      ...STRANDS.map(strand =>
        computeTextWidthInPixels(
          this.config.strandTerminusModifications[strand][STRAND_TO_END_TERMINUS_MAP[strand][end]],
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

  private computeXPositionOfTerminusModifications(): Record<typeof STRAND_ENDS[number], StrandToNumberMap> {
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
    ) as Record<typeof STRAND_ENDS[number], StrandToNumberMap>;

    return xPositionOfTerminusModifications;
  }

  getXPositionOfTerminusModifications(): Record<typeof STRAND_ENDS[number], StrandToNumberMap> {
    return this.xPositionOfTerminusModifications;
  }

  getXPositionOfStrandLabels(): StrandEndToNumberMap {
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
    return computeTotalSVGHeight(this.config.isAntisenseStrandIncluded);
  }

  computeNumericLabelXOffset(bases: string[], nucleobaseIndex: number, nucleotideNumericLabel: number): number {
    const criterion = isSingleDigitNumber(nucleotideNumericLabel) || NUCLEOTIDES.includes(bases[nucleobaseIndex]);
    const shiftAmount = criterion
      ? NUMERIC_LABEL_POSITION_OFFSET.ONE_DIGIT
      : NUMERIC_LABEL_POSITION_OFFSET.TWO_DIGIT;
    return shiftAmount;
  }

  getCenterPositionOfLinkageStar(index: number, strand: STRAND): Position {
    const centerPosition = {
      x: this.getXPositionOfLinkageStar(index, strand),
      y: this.getYPositionOfLinkageStar(strand),
    };
    return centerPosition;
  }

  // todo: remove legacy division of x-y coordinates throughout the code
  private getXPositionOfLinkageStar(index: number, strand: STRAND): number {
    return this._computeNucleobaseCircleXPosition(index, this.rightOverhangs[strand]) + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
  }

  private getYPositionOfLinkageStar(strand: STRAND): number {
    return Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL + SVG_CIRCLE_SIZES.LINKAGE_STAR_RADIUS;
  }

  getCommentLabelPosition(): Position {
    const xPositionOfStrandLabels = this.getXPositionOfStrandLabels();
    const commentLabelPosition = {
      x: xPositionOfStrandLabels[STRAND_END.LEFT],
      y: computeCommentYPosition(this.config.isAntisenseStrandIncluded),
    };
    return commentLabelPosition;
  }

  getStarElementLabelPosition(): Position {
    const starElementLabelPosition = {
      x: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
      y: computeLegendCircleYPosition(this.config.isAntisenseStrandIncluded),
    };

    return starElementLabelPosition;
  }

  getPhosphorothioateLinkageLabelPosition(): Position {
    const position = {
      x: 2 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS - 8,
      y: computeLegendTextYPosition(this.config.isAntisenseStrandIncluded),
    };

    return position;
  }

  getTerminalModificationLabelPosition(strand: STRAND, end: STRAND_END): Position {
    const xPosition = this.getXPositionOfTerminusModifications()[end][strand];
    const labelPosition = {
      x: this.xPositionOfTerminusModifications[end][strand],
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL,
    };

    return labelPosition;
  }

  getStrandEndLabelPosition(strand: STRAND, end: STRAND_END): Position {
    const xPosition = this.getXPositionOfStrandLabels()[end];
    const labelPosition = {
      x: xPosition,
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL,
    };

    return labelPosition;
  }
}
