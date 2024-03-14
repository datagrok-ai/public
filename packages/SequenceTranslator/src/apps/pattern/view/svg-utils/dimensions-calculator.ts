import {NUCLEOTIDES} from '../../../common/model/const';
import { STRAND, STRANDS, STRAND_END, STRAND_TO_END_TERMINUS_MAP, STRAND_ENDS } from '../../model/const';
import {PatternConfiguration} from '../../model/types';
import { SVG_CIRCLE_SIZES, SVG_TEXT_FONT_SIZES, NUMERIC_LABEL_POSITION_OFFSET, DEFAULT_FONT_FAMILY, Y_POSITIONS_FOR_STRAND_ELEMENTS, STRAND_END_LABEL_TEXT} from './const';
import {Position, StrandToNumberMap, StrandEndToNumberMap} from '../types';
import {isOverhangNucleotide} from '../../model/utils';

export class PatternSVGDimensionsCalculator {
  private canvasDimensions: CanvasDimensionCalculator;
  private nucleotidePositionCalculator: NucleotidePositionCalculator;
  private legendPositionCalculator: LegendPositionCalculator;
  private labelPositionCalculator: LabelPositionCalculator;
  private linkageStarPositionCalculator: LinkageStarPositionCalculator;

  constructor(
    private config: PatternConfiguration,
  ) {
    const rightOverhangNucleotideCounts = this.computeRightOverhangNucleotideCounts();
    const maxEffectiveStrandLength = this.computeMaxEffectiveStrandLength(rightOverhangNucleotideCounts);
    const maxWidthOfRightOverhangs = this.computeMaxWidthOfRightOverhangs(rightOverhangNucleotideCounts);
    const strandLabelWidth = this.computeStrandLabelWidth();
    const maxTerminusLabelWidthByEnd = this.computeMaxWidthOfTerminusLabels();

    this.initializeCalculators(
      rightOverhangNucleotideCounts,
      maxEffectiveStrandLength,
      maxWidthOfRightOverhangs,
      strandLabelWidth,
      maxTerminusLabelWidthByEnd,
    );
  }

  private initializeCalculators(
    rightOverhangNucleotideCounts: StrandToNumberMap,
    maxEffectiveStrandLength: number,
    maxWidthOfRightOverhangs: number,
    strandLabelWidth: StrandEndToNumberMap,
    maxTerminusLabelWidthByEnd: StrandEndToNumberMap,
  ): void {
    this.canvasDimensions = new CanvasDimensionCalculator(
      this.config,
      maxEffectiveStrandLength,
      maxWidthOfRightOverhangs,
      strandLabelWidth,
      maxTerminusLabelWidthByEnd,
    );

    this.nucleotidePositionCalculator = new NucleotidePositionCalculator(
      this.config,
      maxEffectiveStrandLength,
      rightOverhangNucleotideCounts,
      maxTerminusLabelWidthByEnd,
    );

    this.legendPositionCalculator = new LegendPositionCalculator(
      this.config,
      this.canvasDimensions,
    );

    this.labelPositionCalculator = new LabelPositionCalculator(
      this.config,
      maxWidthOfRightOverhangs,
      maxTerminusLabelWidthByEnd,
      strandLabelWidth,
      rightOverhangNucleotideCounts,
      this.nucleotidePositionCalculator,
    );

    this.linkageStarPositionCalculator = new LinkageStarPositionCalculator(
      this.nucleotidePositionCalculator,
      rightOverhangNucleotideCounts,
    );
  }

  private computeRightOverhangNucleotideCounts(): StrandToNumberMap {
    return STRANDS.reduce((overhangCounts, strand) => {
      overhangCounts[strand] = this.countOverhangNucleotidesAtStartOfStrand(strand);
      return overhangCounts;
    }, {} as StrandToNumberMap);
  }

  private computeMaxEffectiveStrandLength(rightOverhangNucleotideCounts: StrandToNumberMap): number {
    return Math.max(...STRANDS.map(strand => this.config.nucleotideSequences[strand].length - rightOverhangNucleotideCounts[strand]));
  }

  private computeMaxWidthOfRightOverhangs(rightOverhangNucleotideCounts: StrandToNumberMap): number {
    return Math.max(rightOverhangNucleotideCounts[STRAND.SENSE], rightOverhangNucleotideCounts[STRAND.ANTISENSE]);
  }

  getCanvasWidth(): number {
    return this.canvasDimensions.getCanvasWidth();
  }

  getCanvasHeight(): number {
    return this.canvasDimensions.getCanvasHeight();
  }

  getNucleotideCirclePosition(indexOfNucleotide: number, strand: STRAND): Position {
    return this.nucleotidePositionCalculator.getNucleotideCirclePosition(indexOfNucleotide, strand);
  };

  getNucleotideLabelTextPosition(indexOfNucleotide: number, strand: STRAND): Position {
    return this.nucleotidePositionCalculator.getNucleotideLabelTextPosition(indexOfNucleotide, strand);
  }

  getNumericLabelPosition(
    indexOfNucleotide: number,
    strand: STRAND,
    displayedNucleotideNumber: number
  ): Position {
    return this.nucleotidePositionCalculator.getNumericLabelPosition(indexOfNucleotide, strand, displayedNucleotideNumber);
  }


  getLegendCirclePosition(index: number, distinctNucleobaseTypes: string[], containsPhosphorothioateLinkages: boolean): Position {
    return this.legendPositionCalculator.getLegendCirclePosition(index, distinctNucleobaseTypes, containsPhosphorothioateLinkages);
  }

  getLegendTextPosition(index: number, distinctNucleobaseTypes: string[], containsPhosphorothioateLinkages: boolean): Position {
    return this.legendPositionCalculator.getLegendTextPosition(index, distinctNucleobaseTypes, containsPhosphorothioateLinkages);
  }

  getStarLabelPosition(): Position {
    return this.legendPositionCalculator.getStarLabelPosition();
  }

  getPhosphorothioateLinkageLabelPosition(): Position {
    return this.legendPositionCalculator.getPhosphorothioateLinkageLabelPosition();
  }

  getTerminusLabelPosition(strand: STRAND, end: STRAND_END): Position {
    return this.labelPositionCalculator.getTerminusLabelPosition(strand, end);
  }

  getStrandEndLabelPosition(strand: STRAND, end: STRAND_END): Position {
    return this.labelPositionCalculator.getStrandEndLabelPosition(strand, end);
  }

  getTitleTextPosition(): Position {
    return this.labelPositionCalculator.getTitleTextPosition();
  }

  getCommentLabelPosition(): Position {
    return this.labelPositionCalculator.getCommentLabelPosition();
  }

  getCenterPositionOfLinkageStar(index: number, strand: STRAND): Position {
    return this.linkageStarPositionCalculator.getCenterPositionOfLinkageStar(index, strand);
  }

  private computeStrandLabelWidth(): StrandEndToNumberMap {
    const widthOfStrandLabel = Object.fromEntries(
      STRAND_ENDS.map(
        strandEnd => [strandEnd, this.getMaxWidthStrandEndLabelsByEnd(strandEnd)]
      )
    ) as StrandEndToNumberMap;

    return widthOfStrandLabel;
  }

  private computeMaxWidthOfTerminusLabels(): StrandEndToNumberMap {
    const maxWidthOfTerminusLabelsByEnd = STRAND_ENDS.reduce((maxWidthMap, end) => {
      maxWidthMap[end] = this.getMaxWidthOfTerminusLabelsByEnd(end);
      return maxWidthMap;
    }, {} as StrandEndToNumberMap);

    return maxWidthOfTerminusLabelsByEnd;
  }

  private getMaxWidthOfTerminusLabelsByEnd(end: STRAND_END): number {
    return this.calculateMaxWidthOfStrandEndLabel(
      strand => this.config.strandTerminusModifications[strand][STRAND_TO_END_TERMINUS_MAP[strand][end]]
    );
  }

  private getMaxWidthStrandEndLabelsByEnd(strandEnd: STRAND_END): number {
    return this.calculateMaxWidthOfStrandEndLabel(
      strand => STRAND_END_LABEL_TEXT[strandEnd][strand]
    );
  }

  private calculateMaxWidthOfStrandEndLabel(getLabelText: (strand: STRAND) => string): number {
    const textWidthCalculator = TextWidthCalculator.getInstance();
    return Math.max(
      ...STRANDS.map(strand =>
        textWidthCalculator.computeTextWidth(
          getLabelText(strand),
          SVG_TEXT_FONT_SIZES.NUCLEOBASE,
          DEFAULT_FONT_FAMILY
        )
      )
    );
  }

  private countOverhangNucleotidesAtStartOfStrand(strand: STRAND): number {
    const nucleotides = this.config.nucleotideSequences[strand];

    let overhangNucleotidesCount = 0;
    for (const nucleotide of nucleotides) {
      if (!isOverhangNucleotide(nucleotide))
        break;
      overhangNucleotidesCount++;
    }

    return overhangNucleotidesCount;
  }
}

class TextWidthCalculator {
  private static instance: TextWidthCalculator;
  private canvas: HTMLCanvasElement;
  private context: CanvasRenderingContext2D | null;
  private pixelRatio: number;

  // WARNING: singleton used to avoid creating canvas element on every call
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

class CanvasDimensionCalculator {
  constructor(
    private config: PatternConfiguration,
    private maxEffectiveStrandLength: number,
    private maxWidthOfRightOverhangs: number,
    private strandLabelWidth: StrandEndToNumberMap,
    private maxTerminusWidthByEnd: StrandEndToNumberMap,
  ) {}

  getCanvasWidth(): number {
    const widthOfNucleobases = SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER * (this.maxEffectiveStrandLength + this.maxWidthOfRightOverhangs);

    const canvasWidth = STRAND_ENDS.reduce((totalWidth, end) => {
      totalWidth += this.strandLabelWidth[end] + this.maxTerminusWidthByEnd[end];
      return totalWidth;
    }, 0) + widthOfNucleobases + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;

    return canvasWidth;
  }

  getCanvasHeight(): number {
    return (this.config.isAntisenseStrandIncluded ? 11 : 9) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
  }
}

class NucleotidePositionCalculator {
  constructor(
    private config: PatternConfiguration,
    private maxEffectiveStrandLength: number,
    private rightOverhangNucleotideCounts: StrandToNumberMap,
    private maxTerminusWidthByEnd: StrandEndToNumberMap,
  ) {}

  getNucleotideCirclePosition(indexOfNucleotide: number, strand: STRAND): Position {
    return {
      x: this.computeNucleotideCircleXPositionByStrand(indexOfNucleotide, strand),
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_CIRCLE,
    };
  };

  getNucleotideLabelTextPosition(indexOfNucleotide: number, strand: STRAND): Position {
    return {
      // todo: refactor legacy dependency on one-digit parameter
      x: this.computeNucleotideCircleXPositionByStrand(indexOfNucleotide, strand) + NUMERIC_LABEL_POSITION_OFFSET.ONE_DIGIT,
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL,
    };
  }

  // todo: cleanup legacy logic
  private computeNucleotideCircleXPositionByStrand(index: number, strand: STRAND): number {
    const rightOverhangCount = this.rightOverhangNucleotideCounts[strand];
    return this.computeNucleobaseCircleXPosition(index, rightOverhangCount);
  }

  computeNucleobaseCircleXPosition(index: number, rightOverhangCount: number): number {
    const rightModificationOffset = this.maxTerminusWidthByEnd[STRAND_END.RIGHT];
    const positionalIndex = this.maxEffectiveStrandLength - index + rightOverhangCount + 1;
    const xPosition = positionalIndex * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;
    const finalPosition = rightModificationOffset + xPosition;
    return finalPosition;
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

    const numericLabelXPosition = this.computeNucleotideCircleXPositionByStrand(indexOfNucleotide, strand) + numericLabelOffset;

    return {
      x: numericLabelXPosition,
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUMERIC_LABEL,
    };
  }

  private computeNumericLabelXOffset(bases: string[], nucleobaseIndex: number, nucleotideNumericLabel: number): number {
    const isSingleDigitLabel = nucleotideNumericLabel >= 0 && nucleotideNumericLabel < 10;

    const criterion = isSingleDigitLabel || NUCLEOTIDES.includes(bases[nucleobaseIndex]);

    return criterion ?
      NUMERIC_LABEL_POSITION_OFFSET.ONE_DIGIT :
      NUMERIC_LABEL_POSITION_OFFSET.TWO_DIGIT;
  }

  /** Inverses index for antisense strand image */
  private getVisualStrandIndex(indexOfNucleotide: number, strand: STRAND): number {
    return strand === STRAND.SENSE ?
      indexOfNucleotide :
      this.config.nucleotideSequences[strand].length - indexOfNucleotide;
  }
}

class LegendPositionCalculator {
  constructor(
    private config: PatternConfiguration,
    private canvasDimensionCalculator: CanvasDimensionCalculator,
  ) {}

  getLegendCirclePosition(index: number, distinctNucleobaseTypes: string[], containsPhosphorothioateLinkages: boolean): Position {
    const centerPosition = {
      x: this.computeLegendCircleXPosition(index, distinctNucleobaseTypes, containsPhosphorothioateLinkages),
      y: this.getLegendVerticalPosition(),
    };
    return centerPosition;
  }

  getLegendTextPosition(index: number, distinctNucleobaseTypes: string[], containsPhosphorothioateLinkages: boolean): Position {
    const legendPosition = {
      x: this.computeLegendCircleXPosition(index, distinctNucleobaseTypes, containsPhosphorothioateLinkages) + SVG_CIRCLE_SIZES.LEGEND_RADIUS + 4,
      y: this.getLegendTextVerticalPosition(),
    };
    return legendPosition;
  }

  private computeLegendCircleXPosition(index: number, distinctNucleobaseTypes: string[], containsPhosphorothioateLinkages: boolean): number {
    const legendStartIndex = containsPhosphorothioateLinkages ? 1 : 0;
    const totalPositions = distinctNucleobaseTypes.length + legendStartIndex;
    const width = this.canvasDimensionCalculator.getCanvasWidth();
    const spacingUnit = width / totalPositions;
    const position = (index + legendStartIndex) * spacingUnit;
    const adjustedPosition = position + SVG_CIRCLE_SIZES.LEGEND_RADIUS;
    return Math.round(adjustedPosition);
  }

  getStarLabelPosition(): Position {
    return {
      x: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
      y: this.getLegendVerticalPosition(),
    };
  }

  getPhosphorothioateLinkageLabelPosition(): Position {
    return {
      x: 2 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS - 8,
      y: this.getLegendTextVerticalPosition(),
    };
  }

  private getLegendVerticalPosition(): number {
    return (this.config.isAntisenseStrandIncluded ? 9.5 : 6) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
  }

  private getLegendTextVerticalPosition(): number {
    const position = this.config.isAntisenseStrandIncluded ? 10 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS : Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].NUCLEOBASE_CIRCLE;
    return position - 3;
  }
}

class LabelPositionCalculator {
  private xPositionOfTerminusModifications: Record<typeof STRAND_ENDS[number], StrandToNumberMap>;
  constructor(
    private config: PatternConfiguration,
    private maxWidthOfRightOverhangs: number,
    private maxTerminusWidthByEnd: StrandEndToNumberMap,
    private strandLabelWidth: StrandEndToNumberMap,
    private rightOverhangNucleotideCounts: StrandToNumberMap,
    private nucleotidePositionCalculator: NucleotidePositionCalculator,
  ) {
    this.xPositionOfTerminusModifications = this.computeXPositionOfTerminusModifications();
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
      // todo: remove legacy grouping by y positions
      y: Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL,
    };
  }

  getTitleTextPosition(): Position {
    return {
      x: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
      y: SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
    };
  }

  getCommentLabelPosition(): Position {
    const y = (this.config.isAntisenseStrandIncluded ? 11 : 8.5) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;

    const x = this.getXPositionOfStrandLabels()[STRAND_END.LEFT];
    return { x, y };
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

  private computeXPositionOfTerminusModifications(): Record<typeof STRAND_ENDS[number], StrandToNumberMap> {
    const xPositionOfTerminusModifications = Object.fromEntries(
      STRAND_ENDS.map(end => [
        end,
        Object.fromEntries(
          STRANDS.map(strand => [
            strand,
            end === STRAND_END.LEFT
            ? this.strandLabelWidth[STRAND_END.LEFT] - 5
            : this.rightOverhangNucleotideCounts[strand] * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER + this.nucleotidePositionCalculator.computeNucleobaseCircleXPosition(-0.5, 0)
          ])
        )
      ])
    ) as Record<typeof STRAND_ENDS[number], StrandToNumberMap>;

    return xPositionOfTerminusModifications;
  }

}

class LinkageStarPositionCalculator {
  constructor(
    private nucleotidePositionCalculator: NucleotidePositionCalculator,
    private rightOverhangNucleotideCounts: StrandToNumberMap,
  ) {}

  getCenterPositionOfLinkageStar(index: number, strand: STRAND): Position {
    return {
      x: this.getXPositionOfLinkageStar(index, strand),
      y: this.getYPositionOfLinkageStar(strand),
    };
  }

  // todo: remove legacy division of x-y coordinates throughout the code
  private getXPositionOfLinkageStar(index: number, strand: STRAND): number {
    return this.nucleotidePositionCalculator.computeNucleobaseCircleXPosition(index, this.rightOverhangNucleotideCounts[strand]) + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
  }

  private getYPositionOfLinkageStar(strand: STRAND): number {
    return Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL + SVG_CIRCLE_SIZES.LINKAGE_STAR_RADIUS;
  }
}
