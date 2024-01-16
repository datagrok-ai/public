import {NUCLEOTIDES} from '../../../common/model/const';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../../model/helpers';
import { STRAND, STRANDS, STRAND_END, STRAND_TO_END_TERMINUS_MAP, STRAND_ENDS, TERMINUS, TERMINI } from '../../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../../model/types';
import { SVG_CIRCLE_SIZES, SVG_TEXT_FONT_SIZES, SVG_ELEMENT_COLORS, STRAND_END_LABEL_TEXT, NUMERIC_LABEL_POSITION_OFFSET, DEFAULT_FONT_FAMILY, Y_POSITIONS_FOR_STRAND_ELEMENTS} from './const';
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
import {PatternSVGDimensionsCalculator} from './dimensions-calculator';

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
    const patternSVGCanvas = this.createCanvas();

    const labelElements = this.createLabelElements();
    // todo: remove repeating code
    patternSVGCanvas.append(...labelElements);

    const countOfNucleotidesExcludingOverhangs = this.getCountOfNucleotidesExcludingOverhangs();
    const strandElements = this.createStrandElements(countOfNucleotidesExcludingOverhangs);
    patternSVGCanvas.append(...strandElements);

    const svgDimensionsManager = this.svgDimensionsManager;
    const svgElementFactory = this.svgElementFactory;

    const titleElement = this.createTitleElement(this.config.patternName, countOfNucleotidesExcludingOverhangs, this.config.isAntisenseStrandIncluded);

    patternSVGCanvas.append(titleElement);

    const legend = new LegendManager(this.config, svgDimensionsManager, svgElementFactory).getLegendItems(this.distinctNucleobaseTypes);
    patternSVGCanvas.append(...legend);

    return patternSVGCanvas;
  }

  private createTitleElement(patternName: string,
    numberOfNucleotides: StrandToNumberMap,
    isAntisenseStrandActive: boolean
  ) {
    const titleText = `${patternName} for ${numberOfNucleotides[STRAND.SENSE]}${isAntisenseStrandActive ? `/${numberOfNucleotides[STRAND.ANTISENSE]}` : ''}mer`;

    const titleTextPosition = this.svgDimensionsManager.getTitleTextPosition();
    return this.svgElementFactory.createTextElement(
      titleText,
      titleTextPosition,
      SVG_TEXT_FONT_SIZES.NUCLEOBASE,
      SVG_ELEMENT_COLORS.TITLE_TEXT
    );
  }

  private createPhosphorothioateLinkageStar(strand: StrandType, index: number): SVGElement | null {
    const isActive = this.config.phosphorothioateLinkageFlags[strand][index];
    if (!isActive)
      return null;

    const centerPosition = this.svgDimensionsManager.getCenterPositionOfLinkageStar(index, strand);
    const color = SVG_ELEMENT_COLORS.LINKAGE_STAR;
    const starElement = this.svgElementFactory.createStarElement(centerPosition, color);

    return starElement;
  }

  // todo reduce the # of args
  // port everything that can be ported to external methods
  private createSvgElementsForSpecificNucleotide(
    indexOfNucleotide: number,
    strand: STRAND,
    counter: NucleotideNumericLabelCounter,
    countOfNucleotidesExcludingOverhangs: StrandToNumberMap
  ): SVGElement[] {
    const nucleotide = this.config.nucleotideSequences[strand][indexOfNucleotide];
    const isOverhang = isOverhangNucleotide(nucleotide);

    counter.decrementIfNotOverhang(isOverhang, strand);
    const displayedNucleotideNumber = strand === STRAND.SENSE ?
      counter.getCurrentCount(strand) + 1 :
      countOfNucleotidesExcludingOverhangs[strand] - counter.getCurrentCount(strand);


    const nucleotideNumericLabelElement = this.createNucleotideNumericLabelElement(indexOfNucleotide, strand, displayedNucleotideNumber);
    const nucleotideCircleElement = this.createNucleobaseCircleElement(indexOfNucleotide, strand);
    const nucleotideLabelTextElement = this.createNucleobaseLabelTextElement(indexOfNucleotide, strand);

    const phosphorothioateLinkageStarElement = this.createPhosphorothioateLinkageStar(strand, indexOfNucleotide);

    const lastNucleotideIndex = this.config.nucleotideSequences[strand].length;
    const lastStar = this.createPhosphorothioateLinkageStar(strand, lastNucleotideIndex);

    const nucleotideSvgElements = [
      nucleotideNumericLabelElement,
      nucleotideCircleElement,
      nucleotideLabelTextElement,
      phosphorothioateLinkageStarElement,
      lastStar,
    ].filter((element) => element !== null) as SVGElement[];

    return nucleotideSvgElements;
  }

  private createNucleotideNumericLabelElement(
    indexOfNucleotide: number,
    strand: STRAND,
    displayedNucleotideNumber: number
  ): SVGElement {
    const nucleotide = this.config.nucleotideSequences[strand][indexOfNucleotide];
    const isOverhang = isOverhangNucleotide(nucleotide);
    const numericLabelPosition = this.svgDimensionsManager.getNumericLabelPosition(indexOfNucleotide, strand, displayedNucleotideNumber);

    const nucleotideNumberLabel = (!isOverhang && this.config.nucleotidesWithNumericLabels.includes(nucleotide))
      ? String(displayedNucleotideNumber)
      : '';

    const nucleotideNumericLabelElement = this.svgElementFactory.createTextElement(
      nucleotideNumberLabel, numericLabelPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT
    );

    return nucleotideNumericLabelElement;
  }

  private createNucleobaseCircleElement(indexOfNucleotide: number, strand: STRAND): SVGCircleElement {
    const nucleotide = this.config.nucleotideSequences[strand][indexOfNucleotide];
    const nucleotideCirclePosition = this.svgDimensionsManager.getNucleotideCirclePosition(indexOfNucleotide, strand);
    const nucleotideCircleElement = this.svgElementFactory.createCircleElement(
      nucleotideCirclePosition, SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS, getNucleobaseColorFromStyleMap(nucleotide)
    );
    return nucleotideCircleElement;
  }

  private createNucleobaseLabelTextElement(indexOfNucleotide: number, strand: STRAND): SVGTextElement {
    const nucleotide = this.config.nucleotideSequences[strand][indexOfNucleotide];
    const nucleobaseLabelTextPosition = this.svgDimensionsManager.getNucleotideLabelTextPosition(indexOfNucleotide, strand);
    const nucleotideLabelTextElement = this.svgElementFactory.createTextElement(
      getNucleobaseLabelForCircle(nucleotide), nucleobaseLabelTextPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, computeTextColorForNucleobaseLabel(nucleotide)
    );
    return nucleotideLabelTextElement;
  }

  private getCountOfNucleotidesExcludingOverhangs(): StrandToNumberMap {
    return STRANDS.reduce((acc, strand) => {
      acc[strand] = this.config.nucleotideSequences[strand].filter(value => !isOverhangNucleotide(value)).length;
      return acc;
    }, {} as StrandToNumberMap);
  }

  private createStrandElements(countOfNucleotidesExcludingOverhangs: StrandToNumberMap): SVGElement[] {
    const svgElements = [] as SVGElement[];

    const nucleotideCounter = new NucleotideNumericLabelCounter(countOfNucleotidesExcludingOverhangs);

    STRANDS.forEach(strand => {
      const criterion = strand === STRAND.SENSE || (strand === STRAND.ANTISENSE && this.config.isAntisenseStrandIncluded);
      if (!criterion)
        return;

      this.config.nucleotideSequences[strand].forEach((_, index) => {
        const elements = this.createSvgElementsForSpecificNucleotide(index, strand, nucleotideCounter, countOfNucleotidesExcludingOverhangs);
        svgElements.push(...elements);
      });
    });

    return svgElements;
  }

  private createLabelElements(): SVGElement[] {
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

    return labelElements;
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

  private createStarElementLabel(isPhosphorothioateLinkageActive: boolean): SVGElement | null {
    const starElementLabelPosition = this.svgDimensionsManager.getStarElementLabelPosition();
    const starElementLabel = isPhosphorothioateLinkageActive
      ? this.svgElementFactory.createStarElement(starElementLabelPosition, SVG_ELEMENT_COLORS.LINKAGE_STAR)
      : null;

    return starElementLabel;
  }

  private getPhosphorothioateLinkageLabel(isPhosphorothioateLinkageActive: boolean): SVGElement | null{
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

class LegendManager {
  constructor(
    private config: PatternConfiguration,
    private svgDimensionsManager: PatternSVGDimensionsCalculator,
    private svgElementFactory: SVGElementFactory,
  ) { }

  getLegendItems(distinctNucleobaseTypes: string[]): SVGElement[] {
    const svgElements = [] as SVGElement[];
    distinctNucleobaseTypes.forEach((nucleobaseType, index) => {
      const legendCircle = this.createLegendCircle(nucleobaseType, index);
      const legendText = this.createLegendText(nucleobaseType, index);
      svgElements.push(legendCircle, legendText);
    });
    return svgElements;
  }

  private createLegendCircle(nucleobaseType: string, index: number): SVGCircleElement {
    const radius = SVG_CIRCLE_SIZES.LEGEND_RADIUS;
    const color = getNucleobaseColorFromStyleMap(nucleobaseType);

    const centerPosition = this.svgDimensionsManager.getLegendCirclePosition(index);
    return this.svgElementFactory.createCircleElement(centerPosition, radius, color);
  }

  private createLegendText(nucleobaseType: string, index: number): SVGTextElement {
    const legendPosition = this.svgDimensionsManager.getLegendTextPosition(index);
    return this.svgElementFactory.createTextElement(nucleobaseType, legendPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);
  }
}

class NucleotideNumericLabelCounter {
  private counts: StrandToNumberMap;
  constructor(initialCouts: StrandToNumberMap) {
    // WARNING: to ensure immutability, we need to deep copy the object
    this.counts = JSON.parse(JSON.stringify(initialCouts)) as StrandToNumberMap;
  }

  decrementIfNotOverhang(isOverhang: boolean, strand: STRAND): void {
    if (!isOverhang)
      this.counts[strand]--;
  }

  getCurrentCount(strand: STRAND): number {
    return this.counts[strand];
  }
}
