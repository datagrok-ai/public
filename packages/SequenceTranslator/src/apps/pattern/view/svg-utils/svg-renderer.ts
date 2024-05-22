import * as DG from 'datagrok-api/dg';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../../model/utils';
import {STRAND, STRANDS, STRAND_END, STRAND_ENDS, TERMINUS, TERMINI} from '../../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../../model/types';
import {SVG_CIRCLE_SIZES, SVG_TEXT_FONT_SIZES, SVG_ELEMENT_COLORS, STRAND_END_LABEL_TEXT} from './const';
import {StrandToNumberMap, StrandEndToSVGElementsMap, TerminusToSVGElementMap} from '../types';
import {
  getNucleobaseLabelForCircle,
  computeTextColorForNucleobaseLabel,
  getNucleobaseColorFromStyleMap,
} from './utils';
import {PatternSVGDimensionsCalculator} from './dimensions-calculator';
import { EventBus } from '../../model/event-bus';

export class NucleotidePatternSVGRenderer {
  private config: PatternConfiguration;
  private patternDimensionsCalculator: PatternSVGDimensionsCalculator;
  private svgFactory: SVGElementFactoryWrapper;
  private strandElementManager: StrandElementBuilder;
  private legendBuilder: LegendBuilder;

  constructor(patternConfig: PatternConfiguration, private eventBus: EventBus) {
    this.setupPatternConfig(patternConfig);

    // todo: prefer dependency injection for all these properties
    this.patternDimensionsCalculator = new PatternSVGDimensionsCalculator(this.config);

    this.svgFactory = new SVGElementFactoryWrapper(new SVGElementFactory(), this.config, this.patternDimensionsCalculator,
      this.eventBus);

    this.strandElementManager = new StrandElementBuilder(this.svgFactory, this.config);
    this.legendBuilder = new LegendBuilder(this.svgFactory, this.config);
  }

  renderPattern(): SVGElement {
    const labelElements = this.createLabelElements();

    const countOfNucleotidesExcludingOverhangs = this.countNucleotidesExcludingOverhangs();
    const strandElements = this.strandElementManager.createStrandElements(countOfNucleotidesExcludingOverhangs);
    const titleElement = this.svgFactory.createTitleElement(this.config.patternName, countOfNucleotidesExcludingOverhangs, this.config.isAntisenseStrandIncluded);

    const legend = this.legendBuilder.getLegendItems();

    const patternSVGCanvas = this.svgFactory.createCanvas();
    patternSVGCanvas.append(...labelElements, ...strandElements, titleElement, ...legend);

    return patternSVGCanvas;
  }

  private countNucleotidesExcludingOverhangs(): StrandToNumberMap {
    return STRANDS.reduce((acc, strand) => {
      acc[strand] = this.config.nucleotideSequences[strand].filter((value) => !isOverhangNucleotide(value)).length;
      return acc;
    }, {} as StrandToNumberMap);
  }

  private createLabelElements(): SVGElement[] {
    const strandEndLabels = this.getStrandEndLabels();
    const labelsTerminusModification = this.getTerminusModificationLabels();

    const labelElements = [
      ...STRANDS.flatMap((strand) => [
        ...Object.values(strandEndLabels[strand]),
        ...Object.values(labelsTerminusModification[strand])
      ]),
    ].filter((element) => element !== null) as SVGElement[];

    return labelElements;
  }

  private setupPatternConfig(patternConfig: PatternConfiguration): void {
    // WARNING: to ensure immutability, we need to deep copy the config object
    this.config = JSON.parse(JSON.stringify(patternConfig)) as PatternConfiguration;

    this.config.nucleotideSequences[STRAND.SENSE].reverse();
    this.config.phosphorothioateLinkageFlags[STRAND.SENSE].reverse();
  }

  private getStrandEndLabels() {
    const strandEndLabels = STRANDS.reduce((acc, strand) => {
      acc[strand] = STRAND_ENDS.reduce((endAcc, end) => {
        endAcc[end] = this.svgFactory.createLabelForStrandEnd(strand, end);
        return endAcc;
      }, {} as StrandEndToSVGElementsMap);
      return acc;
    }, {} as Record<typeof STRANDS[number], StrandEndToSVGElementsMap>);

    return strandEndLabels;
  }

  private getTerminusModificationLabels() {
    const labelsTerminusModification = STRANDS.reduce((acc, strand) => {
      acc[strand] = TERMINI.reduce((terminiAcc, terminus) => {
        terminiAcc[terminus] = this.svgFactory.createTerminusModificationLabel(strand, terminus);
        return terminiAcc;
      }, {} as TerminusToSVGElementMap);
      return acc;
    }, {} as Record<typeof STRANDS[number], TerminusToSVGElementMap>);

    return labelsTerminusModification;
  }
}

class LegendBuilder {
  private containsPhosphorothioateLinkages: boolean;
  constructor(
    private svgFactory: SVGElementFactoryWrapper,
    private config: PatternConfiguration
  ) {
    this.containsPhosphorothioateLinkages = this.checkAnyPhosphorothioateLinkages();
  }

  getLegendItems(): SVGElement[] {
    const commentLabel = this.svgFactory.createCommentLabel();
    const nucleotideLegendItems = this.createLegendItemsForNucleotideTypes();
    const phosphorothioateLinkageLegendItem = this.createLegendItemForPhosphorothioateLinkage(this.containsPhosphorothioateLinkages);

    return [commentLabel, ...nucleotideLegendItems, ...phosphorothioateLinkageLegendItem];
  }

  private createLegendItemsForNucleotideTypes(): SVGElement[] {
    const distinctNucleobaseTypes = this.extractNucleotideTypes();
    const svgElements = [] as SVGElement[];
    distinctNucleobaseTypes.forEach((nucleobaseType, index) => {
      const legendCircle = this.svgFactory.createLegendCircle(nucleobaseType, index, distinctNucleobaseTypes, this.containsPhosphorothioateLinkages);
      const legendText = this.svgFactory.createLegendText(nucleobaseType, index, distinctNucleobaseTypes, this.containsPhosphorothioateLinkages);
      svgElements.push(legendCircle, legendText);
    });
    return svgElements;
  }

  private createLegendItemForPhosphorothioateLinkage(containsPhosphorothioateLinkages: boolean): SVGElement[] {
    const starLinkageLegendLabel = this.svgFactory.createLinkageStarLegendLabel(containsPhosphorothioateLinkages);
    const phosphorothioateLinkageLabel = this.svgFactory.createPhosphorothioateLinkageLabel(containsPhosphorothioateLinkages);

    return [starLinkageLegendLabel, phosphorothioateLinkageLabel].filter((element) => element !== null) as SVGElement[];
  }

  private extractNucleotideTypes(): string[] {
    const distinctNucleotides = [...new Set(
      this.config.nucleotideSequences[STRAND.SENSE].concat(
        this.config.isAntisenseStrandIncluded ? this.config.nucleotideSequences[STRAND.ANTISENSE] : []
      )
    )];

    return distinctNucleotides;
  }

  private checkAnyPhosphorothioateLinkages(): boolean {
    return [...this.config.phosphorothioateLinkageFlags[STRAND.SENSE],
      ...(this.config.isAntisenseStrandIncluded ? this.config.phosphorothioateLinkageFlags[STRAND.ANTISENSE] : [])
    ].some((linkage) => linkage);
  }
}

class NucleotideCountTracker {
  private nucleotideCounts: StrandToNumberMap;
  constructor(initialNucleotideCounts: StrandToNumberMap) {
    // WARNING: to ensure immutability, we need to deep copy the object
    this.nucleotideCounts = JSON.parse(JSON.stringify(initialNucleotideCounts)) as StrandToNumberMap;
  }

  decrementIfNotOverhang(isOverhang: boolean, strand: STRAND): void {
    if (!isOverhang)
      this.nucleotideCounts[strand]--;
  }

  getCurrentCount(strand: STRAND): number {
    return this.nucleotideCounts[strand];
  }
}

class SVGElementFactoryWrapper {
  constructor(
    private svgElementFactory: SVGElementFactory,
    private config: PatternConfiguration,
    private dimensionsCalculator: PatternSVGDimensionsCalculator,
    private eventBus: EventBus
  ) { }

  createCanvas(): SVGElement {
    const canvasWidth = this.dimensionsCalculator.getCanvasWidth();
    const canvasHeight = this.dimensionsCalculator.getCanvasHeight();
    const svgCanvas = this.svgElementFactory.createCanvas(canvasWidth, canvasHeight);

    return svgCanvas;
  }

  createTitleElement(patternName: string,
    numberOfNucleotides: StrandToNumberMap,
    isAntisenseStrandActive: boolean
  ) {
    const titleText = this.getTitleText(patternName, numberOfNucleotides, isAntisenseStrandActive);
    const titleTextPosition = this.dimensionsCalculator.getTitleTextPosition();
    return this.svgElementFactory.createTextElement(
      titleText,
      titleTextPosition,
      SVG_TEXT_FONT_SIZES.NUCLEOBASE,
      SVG_ELEMENT_COLORS.TITLE_TEXT
    );
  }

  private getTitleText(
    patternName: string,
    numberOfNucleotides: StrandToNumberMap,
    isAntisenseStrandActive: boolean
  ): string {
    const senseStrandLength = `${numberOfNucleotides[STRAND.SENSE]}`;
    const antisenseStrandLength = isAntisenseStrandActive ? `/${numberOfNucleotides[STRAND.ANTISENSE]}` : '';
    const titleText = `${patternName} for ${senseStrandLength}${antisenseStrandLength}-mer`;

    return titleText;
  }

  createPhosphorothioateLinkageStar(strand: StrandType, index: number): SVGElement | null {
    const isActive = this.config.phosphorothioateLinkageFlags[strand][index];
    const centerPosition = this.dimensionsCalculator.getCenterPositionOfLinkageStar(index, strand);
    const color = SVG_ELEMENT_COLORS.LINKAGE_STAR;
    const starElement = this.svgElementFactory.createStarElement(centerPosition, color, !isActive ? '0.0' : '1.0');
    const strandIdx = strand === STRAND.SENSE ?
      this.config.nucleotideSequences[strand].length - index : index;
    starElement.onclick = (e) => {
      e.stopPropagation();
      e.preventDefault();
      DG.Menu.popup()
      .item(!isActive ? 'Add PTO' : 'Remove PTO', () => {
        this.eventBus.setPhosphorothioateLinkageFlag(strand, strandIdx, !isActive);
      })
      .show();
    };
    
    return starElement;
  }

  createNucleotideNumericLabel(
    indexOfNucleotide: number,
    strand: STRAND,
    displayedNucleotideNumber: number
  ): SVGElement {
    const nucleotide = this.config.nucleotideSequences[strand][indexOfNucleotide];
    const isOverhang = isOverhangNucleotide(nucleotide);
    const labelPosition = this.dimensionsCalculator.getNumericLabelPosition(indexOfNucleotide, strand, displayedNucleotideNumber);

    const labelText = (!isOverhang && this.config.nucleotidesWithNumericLabels.includes(nucleotide)) ?
      String(displayedNucleotideNumber) :
      '';

    const nucleotideNumericLabel = this.svgElementFactory.createTextElement(
      labelText, labelPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT
    );

    return nucleotideNumericLabel;
  }

  createNucleotideCircle(indexOfNucleotide: number, strand: STRAND): SVGCircleElement {
    const nucleotide = this.config.nucleotideSequences[strand][indexOfNucleotide];
    const nucleotideCirclePosition = this.dimensionsCalculator.getNucleotideCirclePosition(indexOfNucleotide, strand);
    const nucleotideCircle = this.svgElementFactory.createCircleElement(
      nucleotideCirclePosition, SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS, getNucleobaseColorFromStyleMap(nucleotide)
    );
    const strandIdx = strand === STRAND.SENSE ?
      this.config.nucleotideSequences[strand].length - 1 - indexOfNucleotide : indexOfNucleotide;
    nucleotideCircle.onclick = (e) => {
      e.stopPropagation();
      e.preventDefault();
      DG.Menu.popup()
      .item('Remove', () => {
        this.eventBus.removeNucleotide(strand, strandIdx);
      })
      .item('Add', () => {
        this.eventBus.addNucleotide(strand, strandIdx);
      })
      .item('Edit', () => {})
      .show();
    };
    return nucleotideCircle;
  }

  createNucleotideNameLabel(indexOfNucleotide: number, strand: STRAND): SVGTextElement {
    const nucleotide = this.config.nucleotideSequences[strand][indexOfNucleotide];
    const nucleobaseLabelTextPosition = this.dimensionsCalculator.getNucleotideLabelTextPosition(indexOfNucleotide, strand);
    const nucleotideLabelText = this.svgElementFactory.createTextElement(
      getNucleobaseLabelForCircle(nucleotide), nucleobaseLabelTextPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, computeTextColorForNucleobaseLabel(nucleotide)
    );
    return nucleotideLabelText;
  }

  createCommentLabel(): SVGElement {
    const commentLabelPosition = this.dimensionsCalculator.getCommentLabelPosition();

    const commentLabel = this.svgElementFactory.createTextElement(
      this.config.patternComment,
      commentLabelPosition,
      SVG_TEXT_FONT_SIZES.COMMENT,
      SVG_ELEMENT_COLORS.TEXT
    );

    return commentLabel;
  }

  createLinkageStarLegendLabel(isPhosphorothioateLinkageActive: boolean): SVGElement | null {
    const starLabelPosition = this.dimensionsCalculator.getStarLabelPosition();
    const starLabel = isPhosphorothioateLinkageActive ?
      this.svgElementFactory.createStarElement(starLabelPosition, SVG_ELEMENT_COLORS.LINKAGE_STAR, '1.0') :
      null;

    return starLabel;
  }

  createPhosphorothioateLinkageLabel(isPhosphorothioateLinkageActive: boolean): SVGElement | null {
    const position = this.dimensionsCalculator.getPhosphorothioateLinkageLabelPosition();
    const phosphorothioateLinkageLabel = isPhosphorothioateLinkageActive ?
      this.svgElementFactory.createTextElement(
        'ps linkage',
        position,
        SVG_TEXT_FONT_SIZES.COMMENT,
        SVG_ELEMENT_COLORS.TEXT) :
      null;

    return phosphorothioateLinkageLabel;
  };

  createLabelForStrandEnd(strand: StrandType, end: STRAND_END): SVGElement | null {
    const isLabelActive = (strand === STRAND.SENSE) || (strand === STRAND.ANTISENSE && this.config.isAntisenseStrandIncluded);
    if (!isLabelActive)
      return null;


    const labelText = STRAND_END_LABEL_TEXT[end][strand];
    const labelPosition = this.dimensionsCalculator.getStrandEndLabelPosition(strand, end);

    return this.svgElementFactory.createTextElement(labelText, labelPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.TEXT);
  }

  createTerminusModificationLabel(strand: StrandType, terminus: TerminalType): SVGElement | null {
    if (strand === STRAND.ANTISENSE && !this.config.isAntisenseStrandIncluded)
      return null;


    const end = (strand === STRAND.SENSE && terminus === TERMINUS.FIVE_PRIME) ||
      (strand === STRAND.ANTISENSE && terminus === TERMINUS.THREE_PRIME) ? STRAND_END.LEFT : STRAND_END.RIGHT;

    const labelText = this.config.strandTerminusModifications[strand][terminus];
    const labelPosition = this.dimensionsCalculator.getTerminusLabelPosition(strand, end);
    return this.svgElementFactory.createTextElement(labelText, labelPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.MODIFICATION_TEXT);
  }

  createLegendCircle(nucleobaseType: string, index: number, distinctNucleobaseTypes: string[], containsPhosphorothioateLinkages: boolean): SVGCircleElement {
    const radius = SVG_CIRCLE_SIZES.LEGEND_RADIUS;
    const color = getNucleobaseColorFromStyleMap(nucleobaseType);

    const centerPosition = this.dimensionsCalculator.getLegendCirclePosition(index, distinctNucleobaseTypes, containsPhosphorothioateLinkages);
    return this.svgElementFactory.createCircleElement(centerPosition, radius, color);
  }

  createLegendText(nucleobaseType: string, index: number, distinctNucleobaseTypes: string[], containsPhosphorothioateLinkages: boolean): SVGTextElement {
    const legendPosition = this.dimensionsCalculator.getLegendTextPosition(index, distinctNucleobaseTypes, containsPhosphorothioateLinkages);
    return this.svgElementFactory.createTextElement(nucleobaseType, legendPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);
  }
}

class StrandElementBuilder {
  constructor(
    private svgFactory: SVGElementFactoryWrapper,
    private config: PatternConfiguration
  ) { }

  createStrandElements(countOfNucleotidesExcludingOverhangs: StrandToNumberMap): SVGElement[] {
    const svgElements = [] as SVGElement[];

    const nucleotideCounter = new NucleotideCountTracker(countOfNucleotidesExcludingOverhangs);

    STRANDS.forEach((strand) => {
      const criterion = strand === STRAND.SENSE || (strand === STRAND.ANTISENSE && this.config.isAntisenseStrandIncluded);
      if (!criterion)
        return;

      this.config.nucleotideSequences[strand].forEach((_, index) => {
        const elements = this.createElementsForNucleotide(index, strand, nucleotideCounter, countOfNucleotidesExcludingOverhangs);
        svgElements.push(...elements);
      });
    });

    return svgElements;
  }

  // todo reduce the # of args
  private createElementsForNucleotide(
    indexOfNucleotide: number,
    strand: STRAND,
    counter: NucleotideCountTracker,
    countOfNucleotidesExcludingOverhangs: StrandToNumberMap
  ): SVGElement[] {
    const nucleotide = this.config.nucleotideSequences[strand][indexOfNucleotide];
    const isOverhang = isOverhangNucleotide(nucleotide);

    counter.decrementIfNotOverhang(isOverhang, strand);
    const displayedNucleotideNumber = strand === STRAND.SENSE ?
      counter.getCurrentCount(strand) + 1 :
      countOfNucleotidesExcludingOverhangs[strand] - counter.getCurrentCount(strand);


    const nucleotideNumericLabel = this.svgFactory.createNucleotideNumericLabel(indexOfNucleotide, strand, displayedNucleotideNumber);
    const nucleotideCircle = this.svgFactory.createNucleotideCircle(indexOfNucleotide, strand);
    const nucleotideNameLabel = this.svgFactory. createNucleotideNameLabel(indexOfNucleotide, strand);

    const phosphorothioateLinkageStar = this.svgFactory.createPhosphorothioateLinkageStar(strand, indexOfNucleotide);

    const lastNucleotideIndex = this.config.nucleotideSequences[strand].length;
    const lastStar = this.svgFactory.createPhosphorothioateLinkageStar(strand, lastNucleotideIndex);

    const nucleotideSvgElements = [
      nucleotideNumericLabel,
      nucleotideCircle,
      nucleotideNameLabel,
      phosphorothioateLinkageStar,
      lastStar,
    ].filter((element) => element !== null) as SVGElement[];

    return nucleotideSvgElements;
  }
}
