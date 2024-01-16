import {NUCLEOTIDES} from '../../../common/model/const';
import { STRAND, STRANDS, STRAND_END, STRAND_TO_END_TERMINUS_MAP, STRAND_ENDS } from '../../model/const';
import {PatternConfiguration} from '../../model/types';
import { SVG_CIRCLE_SIZES, SVG_TEXT_FONT_SIZES, NUMERIC_LABEL_POSITION_OFFSET, DEFAULT_FONT_FAMILY, Y_POSITIONS_FOR_STRAND_ELEMENTS, STRAND_END_LABEL_TEXT} from './const';
import {Position, StrandToNumberMap, StrandEndToNumberMap} from '../types';
import {isOverhangNucleotide} from '../../model/helpers';

export class PatternSVGDimensionsCalculator {
  private strandLabelWidth: StrandEndToNumberMap;
  private maxTerminusWidthByEnd: StrandEndToNumberMap;
  private strandRightOverhangCounts: StrandToNumberMap;
  private maxLengthOfStrands: number;
  private maxWidthOfRightOverhangs: number;
  private xPositionOfTerminusModifications: Record<typeof STRAND_ENDS[number], StrandToNumberMap>;

  constructor(
    private config: PatternConfiguration,
    //todo: remove from constructor, used by only one method
    // private distinctNucleobaseTypes: string[],
    private isPhosphorothioateLinkageActive: boolean,
  ) {
    this.strandLabelWidth = this.getStrandLabelWidth();
    this.maxTerminusWidthByEnd = this.getMaxWidthOfTerminusLabels();


    this.strandRightOverhangCounts = STRANDS.reduce((acc, strand) => {
      acc[strand] = this.countOverhangNucleotidesAtStartOfStrand(strand);
      return acc;
    }, {} as StrandToNumberMap);

    this.maxLengthOfStrands = Math.max(...STRANDS.map(strand => this.config.nucleotideSequences[strand].length - this.strandRightOverhangCounts[strand]));

    this.maxWidthOfRightOverhangs = Math.max(this.strandRightOverhangCounts[STRAND.SENSE], this.strandRightOverhangCounts[STRAND.ANTISENSE]);

    this.xPositionOfTerminusModifications = this.computeXPositionOfTerminusModifications();
  }

  // todo: more descriptive name
  getCanvasWidth(): number {
    const widthOfNucleobases = SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER * (this.maxLengthOfStrands + this.maxWidthOfRightOverhangs);

    const canvasWidth = STRAND_ENDS.reduce((acc, end) => {
      acc += this.strandLabelWidth[end] + this.maxTerminusWidthByEnd[end];
      return acc;
    }, 0) + widthOfNucleobases + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;

    return canvasWidth;
  }

  getCanvasHeight(): number {
    return (this.config.isAntisenseStrandIncluded ? 11 : 9) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
  }

  getCenterPositionOfLinkageStar(index: number, strand: STRAND): Position {
    return {
      x: this.getXPositionOfLinkageStar(index, strand),
      y: this.getYPositionOfLinkageStar(strand),
    };
  }

  getCommentLabelPosition(): Position {
    const y = (this.config.isAntisenseStrandIncluded ? 11 : 8.5) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;

    const x = this.getXPositionOfStrandLabels()[STRAND_END.LEFT];
    return { x, y };
  }

  getStarLabelPosition(): Position {
    return {
      x: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
      y: this.getLegendVerticalPosition(),
    };
  }

  private getLegendVerticalPosition(): number {
    return (this.config.isAntisenseStrandIncluded ? 9.5 : 6) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
  }

  private getLegendTextVerticalPosition(): number {
    const position = this.config.isAntisenseStrandIncluded ? 10 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS : Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].NUCLEOBASE_CIRCLE;
    return position - 3;
  }

  getPhosphorothioateLinkageLabelPosition(): Position {
    return {
      x: 2 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS - 8,
      y: this.getLegendTextVerticalPosition(),
    };
  }

  getTerminusLabelPosition(strand: STRAND, end: STRAND_END): Position {
    return {
      x: this.xPositionOfTerminusModifications[end][strand],
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL,
    };
  }

  getStrandEndLabelPosition(strand: STRAND, end: STRAND_END): Position {
    const xPosition = this.getXPositionOfStrandLabels()[end];
    return {
      x: xPosition,
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL,
    };
  }

  getTitleTextPosition(): Position {
    return {
      x: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
      y: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
    };
  }

  getNumericLabelPosition(
    indexOfNucleotide: number,
    strand: STRAND,
    displayedNucleotideNumber: number
  ): Position {
    const indexForVisualStrand = this.getVisualStrandIndex(indexOfNucleotide, strand);

    const numericLabelOffset = this.computeNumericLabelXOffset(
      this.config.nucleotideSequences[strand],
      indexForVisualStrand,
      displayedNucleotideNumber,
    );

    const numericLabelXPosition = this.computeNucleotideCircleXPosition(indexOfNucleotide, strand) + numericLabelOffset;

    return {
      x: numericLabelXPosition,
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUMERIC_LABEL,
    };
  }

  getNucleotideCirclePosition(indexOfNucleotide: number, strand: STRAND): Position {
    return {
      x: this.computeNucleotideCircleXPosition(indexOfNucleotide, strand),
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_CIRCLE,
    };
  };

  getNucleotideLabelTextPosition(indexOfNucleotide: number, strand: STRAND): Position {
    return {
      // todo: refactor legacy dependency on one-digit parameter
      x: this.computeNucleotideCircleXPosition(indexOfNucleotide, strand) + NUMERIC_LABEL_POSITION_OFFSET.ONE_DIGIT,
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL,
    };
  }

  getLegendCirclePosition(index: number, distinctNucleobaseTypes: string[]): Position {
    const centerPosition = {
      x: this.computeLegendCircleXPosition(index, distinctNucleobaseTypes),
      y: this.getLegendVerticalPosition(),
    };
    return centerPosition;
  }

  getLegendTextPosition(index: number, distinctNucleobaseTypes: string[]): Position {
    const legendPosition = {
      x: this.computeLegendCircleXPosition(index, distinctNucleobaseTypes) + SVG_CIRCLE_SIZES.LEGEND_RADIUS + 4,
      y: this.getLegendTextVerticalPosition(),
    };
    return legendPosition;
  }

  private getStrandLabelWidth(): StrandEndToNumberMap {
    const widthOfStrandLabel = Object.fromEntries(
      STRAND_ENDS.map(
        strandEnd => [strandEnd, this.computeMaxWidthForStrandEnd(strandEnd)]
      )
    ) as StrandEndToNumberMap;

    return widthOfStrandLabel;
  }

  private getMaxWidthOfTerminusLabels(): StrandEndToNumberMap {
    const maxWidthOfTerminusLabelsByEnd = STRAND_ENDS.reduce((acc, end) => {
      acc[end] = this.getMaxWidthOfTerminusLabelsByEnd(end);
      return acc;
    }, {} as StrandEndToNumberMap);

    return maxWidthOfTerminusLabelsByEnd;
  }

  private getMaxWidthOfTerminusLabelsByEnd(end: typeof STRAND_ENDS[number]) {
    const textWidthCalculator = TextWidthCalculator.getInstance();
    return Math.max(
      ...STRANDS.map(strand =>
        textWidthCalculator.computeTextWidth(
          this.config.strandTerminusModifications[strand][STRAND_TO_END_TERMINUS_MAP[strand][end]],
          SVG_TEXT_FONT_SIZES.NUCLEOBASE,
          DEFAULT_FONT_FAMILY
        )
      )
    );
  };

  private computeLegendCircleXPosition(index: number, distinctNucleobaseTypes: string[]): number {
    const legendStartIndex = this.isPhosphorothioateLinkageActive ? 1 : 0;
    const totalPositions = distinctNucleobaseTypes.length + legendStartIndex;
    const width = this.getCanvasWidth();
    const spacingUnit = width / totalPositions;
    const position = (index + legendStartIndex) * spacingUnit;
    const adjustedPosition = position + SVG_CIRCLE_SIZES.LEGEND_RADIUS;
    return Math.round(adjustedPosition);
  }

  // todo: cleanup legacy logic
  private computeNucleotideCircleXPosition(index: number, strand: STRAND): number {
    const rightOverhangs = this.strandRightOverhangCounts[strand];
    return this._computeNucleobaseCircleXPosition(index, rightOverhangs);
  }

  private _computeNucleobaseCircleXPosition(index: number, rightOverhangs: number): number {
    const rightModificationOffset = this.maxTerminusWidthByEnd[STRAND_END.RIGHT];
    const positionalIndex = this.maxLengthOfStrands - index + rightOverhangs + 1;
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
            ? this.strandLabelWidth[STRAND_END.LEFT] - 5
            : this.strandRightOverhangCounts[strand] * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER + this._computeNucleobaseCircleXPosition(-0.5, 0)
          ])
        )
      ])
    ) as Record<typeof STRAND_ENDS[number], StrandToNumberMap>;

    return xPositionOfTerminusModifications;
  }

  private getXPositionOfStrandLabels(): StrandEndToNumberMap {
    const maxRightTerminusModificationShift = Math.max(
      ...STRANDS.map(strand => this.xPositionOfTerminusModifications[STRAND_END.RIGHT][strand])
    );

    const rightEndXPosition = maxRightTerminusModificationShift
      + this.maxTerminusWidthByEnd[STRAND_END.LEFT]
      + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER * this.maxWidthOfRightOverhangs;

    return {
      [STRAND_END.LEFT]: 0,
      [STRAND_END.RIGHT]: rightEndXPosition,
    };
  }

  private computeNumericLabelXOffset(bases: string[], nucleobaseIndex: number, nucleotideNumericLabel: number): number {
    const isSingleDigitLabel = nucleotideNumericLabel >= 0 && nucleotideNumericLabel < 10;

    const criterion = isSingleDigitLabel || NUCLEOTIDES.includes(bases[nucleobaseIndex]);

    return criterion ?
      NUMERIC_LABEL_POSITION_OFFSET.ONE_DIGIT :
      NUMERIC_LABEL_POSITION_OFFSET.TWO_DIGIT;
  }

  // todo: remove legacy division of x-y coordinates throughout the code
  private getXPositionOfLinkageStar(index: number, strand: STRAND): number {
    return this._computeNucleobaseCircleXPosition(index, this.strandRightOverhangCounts[strand]) + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
  }

  private getYPositionOfLinkageStar(strand: STRAND): number {
    return Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL + SVG_CIRCLE_SIZES.LINKAGE_STAR_RADIUS;
  }

  /** Inverses index for antisense strand image */
  private getVisualStrandIndex(indexOfNucleotide: number, strand: STRAND): number {
    return strand === STRAND.SENSE ?
      indexOfNucleotide :
      this.config.nucleotideSequences[strand].length - indexOfNucleotide;
  }

  private countOverhangNucleotidesAtStartOfStrand(strand: STRAND): number {
    const nucleotides = this.config.nucleotideSequences[strand];

    let overhangCount = 0;
    for (const nucleotide of nucleotides) {
      if (!isOverhangNucleotide(nucleotide))
        break;
      overhangCount++;
    }

    return overhangCount;
  }

  private computeMaxWidthForStrandEnd(strandEnd: STRAND_END): number {
    const textWidthCalculator = TextWidthCalculator.getInstance();
    return Math.max(
      ...STRANDS.map(strand =>
        textWidthCalculator.computeTextWidth(
          STRAND_END_LABEL_TEXT[strandEnd][strand],
          SVG_TEXT_FONT_SIZES.NUCLEOBASE,
          DEFAULT_FONT_FAMILY
        )
      )
    );
  }
}

export class TextWidthCalculator {
  private static instance: TextWidthCalculator;
  private canvas: HTMLCanvasElement;
  private context: CanvasRenderingContext2D | null;
  private pixelRatio: number;

  private constructor() {
    this.canvas = document.createElement('canvas');
    this.context = this.canvas.getContext('2d');
    this.pixelRatio = window.devicePixelRatio || 1;

    this.canvas.width *= this.pixelRatio;
    this.canvas.height *= this.pixelRatio;
  }

  public static getInstance(): TextWidthCalculator {
    if (!TextWidthCalculator.instance)
      TextWidthCalculator.instance = new TextWidthCalculator();

    return TextWidthCalculator.instance;
  }

  public computeTextWidth(text: string, fontSize: number, fontFamily: string): number {
    if (this.context) {
      this.context.font = `${fontSize * this.pixelRatio}px ${fontFamily}`;
      const metrics = this.context.measureText(text);
      return metrics.width / this.pixelRatio;
    }

    return 0;
  }
}
