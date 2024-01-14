import {NUCLEOTIDES} from '../../common/model/const';
import {AXOLABS_STYLE_MAP as styleMap} from '../../common/data-loader/json-loader';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../model/helpers';
import { STRAND, STRANDS, TERMINUS, TERMINI } from '../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../model/types';

const enum STRAND_END {
  LEFT,
  RIGHT,
};

const STRAND_ENDS = [STRAND_END.LEFT, STRAND_END.RIGHT] as const;

const enum LUMINANCE_COEFFICIENTS {
  RED = 0.299,
  GREEN = 0.587,
  BLUE = 0.114,
  THRESHOLD = 186,
};

const enum TEXT_COLOR {
  DARK = '#333333',
  LIGHT = '#ffffff',
};

const enum SVG_CIRCLE_SIZES {
  NUCLEOBASE_RADIUS = 15,
  NUCLEOBASE_DIAMETER = 2 * NUCLEOBASE_RADIUS,
  LEGEND_RADIUS = 6,
  LINKAGE_STAR_RADIUS = 5,
};

const enum SVG_TEXT_FONT_SIZES {
  NUCLEOBASE = 17,
  COMMENT = 14,
};

const enum SVG_ELEMENT_COLORS {
  LINKAGE_STAR = 'red',
  TEXT = 'var(--grey-6)',
  TITLE_TEXT = 'black',
  MODIFICATION_TEXT = 'red'
};

const STRAND_END_LABEL_TEXT = {
  [STRAND_END.LEFT]: {
    [STRAND.SENSE]: `${STRAND.SENSE}: ${TERMINUS.FIVE_PRIME}`,
    [STRAND.ANTISENSE]: `${STRAND.ANTISENSE}: ${TERMINUS.THREE_PRIME}`,
  },
  [STRAND_END.RIGHT]: {
    [STRAND.SENSE]: `${TERMINUS.THREE_PRIME}`,
    [STRAND.ANTISENSE]: `${TERMINUS.FIVE_PRIME}`,
  }
} as const;

const NUMERIC_LABEL_POSITION_OFFSET = {
  ONE_DIGIT: -5,
  TWO_DIGIT: -10,
} as const;

const DEFAULT_FONT_FAMILY = 'Arial';

const Y_POSITIONS_FOR_STRAND_ELEMENTS = {
  [STRAND.SENSE]: {
    NUMERIC_LABEL: 2 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
    NUCLEOBASE_CIRCLE: 3.5 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
    NUCLEOBASE_LABEL: 4 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
  },
  [STRAND.ANTISENSE]: {
    NUMERIC_LABEL: 8.5 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
    NUCLEOBASE_CIRCLE: 6.5 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
    NUCLEOBASE_LABEL: 7 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS,
  }
};

function computeCommentYPosition(isAntisenseStrandActive: boolean): number {
  return (isAntisenseStrandActive ? 11 : 8.5) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
}

function computeLegendCircleYPosition(isAntisenseStrandActive: boolean): number {
  return (isAntisenseStrandActive ? 9.5 : 6) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
}

function computeLegendTextYPosition(isAntisenseStrandActive: boolean): number {
  const position = isAntisenseStrandActive ? 10 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS : Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].NUCLEOBASE_CIRCLE;
  return position - 3;
}

function computeTotalSVGHeight(isAntisenseStrandActive: boolean): number {
  return (isAntisenseStrandActive ? 11 : 9) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
}

export function renderNucleotidePattern(patternConfiguration: PatternConfiguration): Element {
  let {
    patternName,
    isAntisenseStrandActive,
    // todo: rename
    bases: nucleobases,
    phosphorothioateLinkages,
    terminalModifications,
    comment,
    modificationsWithNumericLabels
  } = patternConfiguration;

  function computeLegendCircleXPosition(index: number): number {
    const totalPositions = distinctNucleobaseTypes.length + legendStartIndex;
    const spacingUnit = width / totalPositions;
    const position = (index + legendStartIndex) * spacingUnit;
    const adjustedPosition = position + SVG_CIRCLE_SIZES.LEGEND_RADIUS;
    return Math.round(adjustedPosition);
  }

  function computeNucleobaseCircleXPosition(index: number, rightOverhangs: number): number {
    const rightModificationOffset = maxWidthOfTerminusLabelsByEnd[STRAND_END.RIGHT];
    const positionalIndex = maxStrandLength - index + rightOverhangs + 1;
    const xPosition = positionalIndex * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;
    const finalPosition = rightModificationOffset + xPosition;
    return finalPosition;
  }

  function computeNumericLabelXOffset(nucleobases: string[], nucleobaseIndex: number, nucleotideNumericLabel: number): number {
    const isSingleDigitOrSimpleNucleotide = isSingleDigitNumber(nucleotideNumericLabel) || NUCLEOTIDES.includes(nucleobases[nucleobaseIndex]);
    const shiftAmount = isSingleDigitOrSimpleNucleotide
      ? NUMERIC_LABEL_POSITION_OFFSET.ONE_DIGIT
      : NUMERIC_LABEL_POSITION_OFFSET.TWO_DIGIT;
    return shiftAmount;
  }

  nucleobases[STRAND.SENSE] = nucleobases[STRAND.SENSE].reverse();
  phosphorothioateLinkages[STRAND.SENSE] = phosphorothioateLinkages[STRAND.SENSE].reverse();

  const rightOverhangs = STRANDS.reduce((acc, strand) => {
    acc[strand] = countOverhangNucleotidesAtStrandEnd(nucleobases[strand]);
    return acc;
  }, {} as Record<typeof STRANDS[number], number>);

  const maxStrandLength = Math.max(...STRANDS.map(strand => nucleobases[strand].length - rightOverhangs[strand]));

  const widthOfRightOverhangs = Math.max(rightOverhangs[STRAND.SENSE], rightOverhangs[STRAND.ANTISENSE]);
  const widthOfBases = SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER * (maxStrandLength + widthOfRightOverhangs);

  const strandEndToTerminusMapping = {
    [STRAND.SENSE]: {
      [STRAND_END.LEFT]: TERMINUS.THREE_PRIME,
      [STRAND_END.RIGHT]: TERMINUS.FIVE_PRIME
    },
    [STRAND.ANTISENSE]: {
      [STRAND_END.LEFT]: TERMINUS.FIVE_PRIME,
      [STRAND_END.RIGHT]: TERMINUS.THREE_PRIME
    }
  } as const;

  function computeMaxWidthOfTerminalLabels(end: typeof STRAND_ENDS[number]) {
    return Math.max(
      ...STRANDS.map(strand =>
        computeTextWidthInPixels(
          terminalModifications[strand][strandEndToTerminusMapping[strand][end]],
          SVG_TEXT_FONT_SIZES.NUCLEOBASE,
          DEFAULT_FONT_FAMILY
        )
      )
    );
  };

  const maxWidthOfTerminusLabelsByEnd = STRAND_ENDS.reduce((acc, end) => {
    acc[end] = computeMaxWidthOfTerminalLabels(end);
    return acc;
  }, {} as Record<typeof STRAND_ENDS[number], number>);

  const distinctNucleobaseTypes = [...new Set(
    nucleobases[STRAND.SENSE].concat(
      isAntisenseStrandActive ? nucleobases[STRAND.ANTISENSE] : []
    )
  )];

  const isPhosphorothioateLinkageActive = [
    ...phosphorothioateLinkages[STRAND.SENSE],
    ...(isAntisenseStrandActive ? phosphorothioateLinkages[STRAND.ANTISENSE] : [])
  ].some(linkage => linkage);

  const legendStartIndex = isPhosphorothioateLinkageActive ? 1 : 0;

  const widthOfStrandLabel = Object.fromEntries(
    STRAND_ENDS.map(
      strandEnd => [strandEnd, computeMaxWidthForStrandEnd(strandEnd)]
    )
  ) as Record<typeof STRAND_ENDS[number], number>;

  const xPositionOfTerminusModifications = Object.fromEntries(
    STRAND_ENDS.map(end => [
      end,
      Object.fromEntries(
        STRANDS.map(strand => [
          strand,
          end === STRAND_END.LEFT
            ? widthOfStrandLabel[STRAND_END.LEFT] - 5
            : rightOverhangs[strand] * SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER + computeNucleobaseCircleXPosition(-0.5, 0)
        ])
      )
    ])
  );

  const computeRightEndPosition = (): number => {
    const maxRightTerminusModificationShift = Math.max(
      ...STRANDS.map(strand => xPositionOfTerminusModifications[STRAND_END.RIGHT][strand])
    );
    return maxRightTerminusModificationShift
      + maxWidthOfTerminusLabelsByEnd[STRAND_END.LEFT]
      + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER * widthOfRightOverhangs;
  };

  const xPositionOfStrandLabels = {
    [STRAND_END.LEFT]: 0,
    [STRAND_END.RIGHT]: computeRightEndPosition()
  };

  const width = STRAND_ENDS.reduce((acc, end) => {
    acc += widthOfStrandLabel[end] + maxWidthOfTerminusLabelsByEnd[end];
    return acc;
  }, 0) + widthOfBases + SVG_CIRCLE_SIZES.NUCLEOBASE_DIAMETER;

  const svgElementFactory = new SVGElementFactory();

  const image = svgElementFactory.createCanvas(width, computeTotalSVGHeight(isAntisenseStrandActive));

  function createLabelForStrandEnd(strand: StrandType, end: STRAND_END): SVGElement | null {
    const isLabelActive = (strand === STRAND.SENSE) || (strand === STRAND.ANTISENSE && isAntisenseStrandActive);
    if (!isLabelActive) {
      return null;
    }

    const labelText = STRAND_END_LABEL_TEXT[end][strand];
    const xPosition = xPositionOfStrandLabels[end];
    const yPosition = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL;

    return svgElementFactory.createTextElement(labelText, xPosition, yPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.TEXT);
  }

  const labelsStrandEnd = STRANDS.reduce((acc, strand) => {
    acc[strand] = STRAND_ENDS.reduce((endAcc, end) => {
      endAcc[end] = createLabelForStrandEnd(strand, end);
      return endAcc;
    }, {} as Record<typeof STRAND_ENDS[number], SVGElement | null>);
    return acc;
  }, {} as Record<typeof STRANDS[number], Record<typeof STRAND_ENDS[number], SVGElement | null>>);

  function createTerminusModificationLabel(strand: StrandType, terminus: TerminalType): SVGElement | null {
    if (strand === STRAND.ANTISENSE && !isAntisenseStrandActive) {
      return null;
    }

    const end = (strand === STRAND.SENSE && terminus === TERMINUS.FIVE_PRIME) ||
      (strand === STRAND.ANTISENSE && terminus === TERMINUS.THREE_PRIME) ? STRAND_END.LEFT : STRAND_END.RIGHT;

    const labelText = terminalModifications[strand][terminus];
    const xPosition = xPositionOfTerminusModifications[end][strand];
    const yPosition = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL;

    return svgElementFactory.createTextElement(labelText, xPosition, yPosition, SVG_TEXT_FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.MODIFICATION_TEXT);
  }

  const labelsTerminusModification = STRANDS.reduce((acc, strand) => {
    acc[strand] = TERMINI.reduce((terminiAcc, terminus) => {
      terminiAcc[terminus] = createTerminusModificationLabel(strand, terminus);
      return terminiAcc;
    }, {} as Record<typeof TERMINI[number], SVGElement | null>);
    return acc;
  }, {} as Record<typeof STRANDS[number], Record<typeof TERMINI[number], SVGElement | null>>);

  const commentLabel = svgElementFactory.createTextElement(comment, xPositionOfStrandLabels[STRAND_END.LEFT], computeCommentYPosition(isAntisenseStrandActive), SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);

  const starElement = isPhosphorothioateLinkageActive
    ? svgElementFactory.createStarElement(SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS, computeLegendCircleYPosition(isAntisenseStrandActive), SVG_ELEMENT_COLORS.LINKAGE_STAR)
    : null;

  const psLinkageLabel = isPhosphorothioateLinkageActive ? svgElementFactory.createTextElement('ps linkage', 2 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS - 8, computeLegendTextYPosition(isAntisenseStrandActive), SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT) : null;

  const svgElements = [
    ...STRANDS.map(strand => Object.values(labelsStrandEnd[strand])),
    ...STRANDS.map(strand => Object.values(labelsTerminusModification[strand])),
  ].flat()
    .concat(commentLabel, starElement, psLinkageLabel)
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

    const xPosition = computeNucleobaseCircleXPosition(index, rightOverhangs[strand]) + SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
    const yPosition = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].NUCLEOBASE_LABEL + SVG_CIRCLE_SIZES.LINKAGE_STAR_RADIUS;
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
    const numericLabelOffset = computeNumericLabelXOffset(
      nucleobases[strand],
      nucleobaseIndex,
      displayedNumericLabelIndex,
    );
    const numericLabelXPosition = computeNucleobaseCircleXPosition(index, rightOverhangs[strand]) + numericLabelOffset;

    const nucleotideNumberLabel = (!isOverhang && modificationsWithNumericLabels.includes(nucleobase))
      ? String(displayedNumericLabelIndex)
      : '';

    const nucleotideNumberTextElement = svgElementFactory.createTextElement(
      nucleotideNumberLabel, numericLabelXPosition, yPositions.NUMERIC_LABEL, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT
    );

    const nucleobaseCircleXPosition = computeNucleobaseCircleXPosition(index, rightOverhangs[strand]);
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

  function createLegendCircle(baseType: string, index: number, isAntisenseStrandActive: boolean): SVGCircleElement {
    const xPosition = computeLegendCircleXPosition(index);
    const yPosition = computeLegendCircleYPosition(isAntisenseStrandActive);
    const color = getNucleobaseColorFromStyleMap(baseType);
    return svgElementFactory.createCircleElement(xPosition, yPosition, SVG_CIRCLE_SIZES.LEGEND_RADIUS, color);
  }

  function createLegendText(baseType: string, index: number, isAntisenseStrandActive: boolean): SVGTextElement {
    const xPosition = computeLegendCircleXPosition(index) + SVG_CIRCLE_SIZES.LEGEND_RADIUS + 4;
    const yPosition = computeLegendTextYPosition(isAntisenseStrandActive);
    return svgElementFactory.createTextElement(baseType, xPosition, yPosition, SVG_TEXT_FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);
  }

  distinctNucleobaseTypes.forEach((baseType, index) => {
    const legendCircle = createLegendCircle(baseType, index, isAntisenseStrandActive);
    const legendText = createLegendText(baseType, index, isAntisenseStrandActive);
    image.append(legendCircle, legendText);
  });

  return image;
}

function isSingleDigitNumber(n: number): boolean {
  return n >= 0 && n < 10;
}

function countOverhangNucleotidesAtStrandEnd(modifications: string[]): number {
  const lastIdx = modifications.length - 1;
  let count = 0;
  while (count <= lastIdx && isOverhangNucleotide(modifications[count])) {
    count++;
  }
  return count === lastIdx + 1 ? 0 : count;
}

function computeTextWidthInPixels(text: string, fontSize: number, fontFamily: string): number {
  const canvas = document.createElement('canvas');
  const context = canvas.getContext('2d');
  if (context) {
    context.font = `${fontSize}px ${fontFamily}`;
    const metrics = context.measureText(text);
    return 2 * metrics.width;
  }
  return 0;
}

function getNucleobaseLabelForCircle(nucleobase: string): string {
  const criterion = !isOverhangNucleotide(nucleobase) && NUCLEOTIDES.includes(nucleobase);

  return criterion ? nucleobase : '';
}

function computeTextColorForNucleobaseLabel(nucleobase: string): string {
  const nucleobaseColor = styleMap[nucleobase]?.color || '';

  const rgbValues = nucleobaseColor.match(/\d+/g)?.map(Number);
  if (!rgbValues || rgbValues.length < 3) {
    return TEXT_COLOR.LIGHT;
  }

  const [r, g, b] = rgbValues;
  const luminance = r * LUMINANCE_COEFFICIENTS.RED + g * LUMINANCE_COEFFICIENTS.GREEN + b * LUMINANCE_COEFFICIENTS.BLUE;
  return luminance > LUMINANCE_COEFFICIENTS.THRESHOLD ? TEXT_COLOR.DARK : TEXT_COLOR.LIGHT;
}

function getNucleobaseColorFromStyleMap(nucleobase: string): string {
  return styleMap[nucleobase].color;
}

function computeMaxWidthForStrandEnd(strandEnd: STRAND_END): number {
  return Math.max(
    ...STRANDS.map(strand =>
      computeTextWidthInPixels(
        STRAND_END_LABEL_TEXT[strandEnd][strand],
        SVG_TEXT_FONT_SIZES.NUCLEOBASE,
        DEFAULT_FONT_FAMILY
      )
    )
  );
}
