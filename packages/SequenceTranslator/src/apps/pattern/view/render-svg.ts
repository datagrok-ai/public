import {NUCLEOTIDES} from '../../common/model/const';
import {axolabsStyleMap as styleMap} from '../../common/data-loading-utils/json-loader';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../model/helpers';
import { STRAND, STRANDS, TERMINUS, TERMINI } from '../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../model/types';

const enum STRAND_END {
  LEFT,
  RIGHT,
};

const STRAND_ENDS = [STRAND_END.LEFT, STRAND_END.RIGHT] as const;

const enum LUMINANCE {
  RED_COEFFICIENT = 0.299,
  GREEN_COEFFICIENT = 0.587,
  BLUE_COEFFICIENT = 0.114,
  THRESHOLD = 186,
};

const enum TEXT_COLOR {
  DARK = '#333333',
  LIGHT = '#ffffff',
};

const enum CIRCLE_DIMENSIONS {
  NUCLEOBASE_RADIUS = 15,
  NUCLEOBASE_DIAMETER = 2 * NUCLEOBASE_RADIUS,
  LEGEND_RADIUS = 6,
  LINKAGE_STAR_RADIUS = 5,
};

const enum FONT_SIZES {
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


const TITLE_POSITION = {
  X: CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
  Y: CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
};

const Y_POSITIONS_FOR_STRAND_ELEMENTS = {
  [STRAND.SENSE]: {
    NUMERIC_LABEL: 2 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
    NUCLEOBASE_CIRCLE: 3.5 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
    BASE_LABEL: 4 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
  },
  [STRAND.ANTISENSE]: {
    NUMERIC_LABEL: 8.5 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
    NUCLEOBASE_CIRCLE: 6.5 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
    BASE_LABEL: 7 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
  }
};

const Y_POSITIONS_FOR_SVG_ELEMENTS = {
  COMMENT: (asExists: boolean) => (asExists) ? 11 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS : 8.5 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
  LEGEND_CIRCLES: (asExists: boolean) => (asExists) ? 9.5 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS : 6 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
  LEGEND_TEXT: (asExists: boolean) => (asExists) ? 10 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS - 3 : Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].NUCLEOBASE_CIRCLE - 3,
  SVG_TOTAL_HEIGHT: (asExists: boolean) => (asExists) ? 11 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS : 9 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS,
};

export function renderNucleotidePattern(patternConfiguration: PatternConfiguration): Element {
  let {
    patternName,
    isAntisenseStrandActive,
    bases,
    phosphorothioateLinkages,
    terminalModifications,
    comment,
    modificationsWithNumericLabels
  } = patternConfiguration;

  function computeLegendCircleXPosition(index: number): number {
    const totalPositions = distinctNucleobaseTypes.length + legendStartIndex;
    const spacingUnit = width / totalPositions;
    const position = (index + legendStartIndex) * spacingUnit;
    const adjustedPosition = position + CIRCLE_DIMENSIONS.LEGEND_RADIUS;
    return Math.round(adjustedPosition);
  }

  function computeNucleobaseCircleXPosition(index: number, rightOverhangs: number): number {
    const rightModificationOffset = maxWidthOfTerminusLabelsByEnd[STRAND_END.RIGHT];
    const positionalIndex = maxStrandLength - index + rightOverhangs + 1;
    const xPosition = positionalIndex * CIRCLE_DIMENSIONS.NUCLEOBASE_DIAMETER;
    const finalPosition = rightModificationOffset + xPosition;
    return finalPosition;
  }

  function computeNumericLabelOffset(bases: string[], generalIndex: number, nucleotideIndex: number): number {
    const isSingleDigitOrValidNucleotide = isSingleDigitNumber(nucleotideIndex) || NUCLEOTIDES.includes(bases[generalIndex]);
    const shiftAmount = isSingleDigitOrValidNucleotide
      ? NUMERIC_LABEL_POSITION_OFFSET.ONE_DIGIT
      : NUMERIC_LABEL_POSITION_OFFSET.TWO_DIGIT;
    return shiftAmount;
  }

  bases[STRAND.SENSE] = bases[STRAND.SENSE].reverse();
  phosphorothioateLinkages[STRAND.SENSE] = phosphorothioateLinkages[STRAND.SENSE].reverse();

  const rightOverhangs = STRANDS.reduce((acc, strand) => {
    acc[strand] = countOverhangNucleotidesAtStrandEnd(bases[strand]);
    return acc;
  }, {} as Record<typeof STRANDS[number], number>);

  const maxStrandLength = Math.max(...STRANDS.map(strand => bases[strand].length - rightOverhangs[strand]));

  const widthOfRightOverhangs = Math.max(rightOverhangs[STRAND.SENSE], rightOverhangs[STRAND.ANTISENSE]);
  const widthOfBases = CIRCLE_DIMENSIONS.NUCLEOBASE_DIAMETER * (maxStrandLength + widthOfRightOverhangs);

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
          FONT_SIZES.NUCLEOBASE,
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
    bases[STRAND.SENSE].concat(
      isAntisenseStrandActive ? bases[STRAND.ANTISENSE] : []
    )
  )];

  const isPhosphorothioateLinkageActive = isAntisenseStrandActive
    ? [
        ...phosphorothioateLinkages[STRAND.SENSE],
        ...phosphorothioateLinkages[STRAND.ANTISENSE]
      ].some(linkage => linkage)
    : phosphorothioateLinkages[STRAND.SENSE].some(linkage => linkage);

  const legendStartIndex = isPhosphorothioateLinkageActive ? 1 : 0;

  const widthOfStrandLabel = Object.fromEntries(
    STRAND_ENDS.map(
      strandEnd => [strandEnd, computeMaxWidthForStrandEnd(strandEnd)]
    )
  ) as Record<typeof STRAND_ENDS[number], number>;

  const xPositionOfTerminusModifications = {} as Record<typeof STRAND_ENDS[number], Record<typeof STRANDS[number], number>>;
  xPositionOfTerminusModifications[STRAND_END.LEFT] = STRANDS.reduce((acc, strand) => {
    acc[strand] = widthOfStrandLabel[STRAND_END.LEFT] - 5;
    return acc;
  }, {} as Record<typeof STRANDS[number], number>);
  xPositionOfTerminusModifications[STRAND_END.RIGHT] = STRANDS.reduce((acc, strand) => {
    acc[strand] = rightOverhangs[strand] * CIRCLE_DIMENSIONS.NUCLEOBASE_DIAMETER + computeNucleobaseCircleXPosition(-0.5, 0);
    return acc;
  }, {} as Record<typeof STRANDS[number], number>);

  const xPositionOfStrandLabels = {
    [STRAND_END.LEFT]: 0,
    [STRAND_END.RIGHT]: Math.max(...STRANDS.map(strand => xPositionOfTerminusModifications[STRAND_END.RIGHT][strand]))
      + maxWidthOfTerminusLabelsByEnd[STRAND_END.LEFT]
      + CIRCLE_DIMENSIONS.NUCLEOBASE_DIAMETER * widthOfRightOverhangs
  };

  const width = widthOfStrandLabel[STRAND_END.LEFT] + maxWidthOfTerminusLabelsByEnd[STRAND_END.LEFT] + widthOfBases + maxWidthOfTerminusLabelsByEnd[STRAND_END.RIGHT] + widthOfStrandLabel[STRAND_END.RIGHT] + CIRCLE_DIMENSIONS.NUCLEOBASE_DIAMETER;

  const svgElementFactory = new SVGElementFactory();

  const image = svgElementFactory.createCanvas(width, Y_POSITIONS_FOR_SVG_ELEMENTS.SVG_TOTAL_HEIGHT(isAntisenseStrandActive));

  const labelsStrandEnd = {} as Record<typeof STRANDS[number], Record<typeof STRAND_ENDS[number], SVGElement | null>>;
  STRANDS.forEach((strand) => {
    labelsStrandEnd[strand] = {} as Record<typeof STRAND_ENDS[number], SVGElement | null>;
    STRAND_ENDS.forEach((end) => {
      const isLabelActive = (strand === STRAND.SENSE) || (strand === STRAND.ANTISENSE && isAntisenseStrandActive);
      const labelText = STRAND_END_LABEL_TEXT[end][strand];
      const xPosition = xPositionOfStrandLabels[end];
      const labelY = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].BASE_LABEL;

      labelsStrandEnd[strand][end] = isLabelActive
        ? svgElementFactory.createTextElement(labelText, xPosition, labelY, FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.TEXT)
        : null;
    });
  });

  const labelsTerminusModification = {} as Record<typeof STRANDS[number], Record<typeof TERMINI[number], SVGElement | null>>;

  STRANDS.forEach((strand) => {
    labelsTerminusModification[strand] = {} as Record<typeof TERMINI[number], SVGElement | null>;

    TERMINI.forEach((terminus) => {
      if (strand === STRAND.ANTISENSE && !isAntisenseStrandActive) {
        labelsTerminusModification[strand][terminus] = null;
        return;
      }

      const end = (strand === STRAND.SENSE && terminus === TERMINUS.FIVE_PRIME) ||
        (strand === STRAND.ANTISENSE && terminus === TERMINUS.THREE_PRIME) ? STRAND_END.LEFT : STRAND_END.RIGHT;
      const labelY = Y_POSITIONS_FOR_STRAND_ELEMENTS[strand].BASE_LABEL;
      const text = terminalModifications[strand][terminus];
      const xPosition = xPositionOfTerminusModifications[end][strand];

      labelsTerminusModification[strand][terminus] = svgElementFactory.createTextElement(
        text, xPosition, labelY, FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.MODIFICATION_TEXT
      );
    });
  });

  const commentLabel = svgElementFactory.createTextElement(comment, xPositionOfStrandLabels[STRAND_END.LEFT], Y_POSITIONS_FOR_SVG_ELEMENTS.COMMENT(isAntisenseStrandActive), FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT);

  const starElement = isPhosphorothioateLinkageActive ? svgElementFactory.createStarElement(CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS, Y_POSITIONS_FOR_SVG_ELEMENTS.LEGEND_CIRCLES(isAntisenseStrandActive), SVG_ELEMENT_COLORS.LINKAGE_STAR) : null;

  const psLinkageLabel = isPhosphorothioateLinkageActive ? svgElementFactory.createTextElement('ps linkage', 2 * CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS - 8, Y_POSITIONS_FOR_SVG_ELEMENTS.LEGEND_TEXT(isAntisenseStrandActive), FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT) : null;

  const svgElements = [
    labelsStrandEnd[STRAND.SENSE][STRAND_END.LEFT],
    labelsStrandEnd[STRAND.SENSE][STRAND_END.RIGHT],
    labelsStrandEnd[STRAND.ANTISENSE][STRAND_END.LEFT],
    labelsStrandEnd[STRAND.ANTISENSE][STRAND_END.RIGHT],
    labelsTerminusModification[STRAND.SENSE][TERMINUS.FIVE_PRIME],
    labelsTerminusModification[STRAND.ANTISENSE][TERMINUS.THREE_PRIME],
    labelsTerminusModification[STRAND.SENSE][TERMINUS.THREE_PRIME],
    labelsTerminusModification[STRAND.ANTISENSE][TERMINUS.FIVE_PRIME],
    commentLabel,
    starElement,
    psLinkageLabel,
  ].filter((element) => element !== null) as SVGElement[];

  image.append(...svgElements);

  const numberOfSsNucleotides = bases[STRAND.SENSE].filter((value) => !isOverhangNucleotide(value)).length;
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = bases[STRAND.SENSE].length - 1; i > -1; i--) {
    const xOfNumbers = computeNucleobaseCircleXPosition(i, rightOverhangs[STRAND.SENSE]) +
      computeNumericLabelOffset(bases[STRAND.SENSE], bases[STRAND.SENSE].length - i, numberOfSsNucleotides - nucleotideCounter);
    if (!isOverhangNucleotide(bases[STRAND.SENSE][i]))
      nucleotideCounter--;
    const nucleotideLabel = (!isOverhangNucleotide(bases[STRAND.SENSE][i]) && modificationsWithNumericLabels.includes(bases[STRAND.SENSE][i])) ?
      String(numberOfSsNucleotides - nucleotideCounter) : '';

    const baseCircleXCoord = computeNucleobaseCircleXPosition(i, rightOverhangs[STRAND.SENSE]);

    const nucleotideTextLabel = svgElementFactory.createTextElement(
      nucleotideLabel, xOfNumbers, Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.SENSE].NUMERIC_LABEL, FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT
    );

    const nucleotideCircle = svgElementFactory.createCircleElement(
      baseCircleXCoord, Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.SENSE].NUCLEOBASE_CIRCLE, CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS, retrieveNucleobaseColor(bases[STRAND.SENSE][i])
    );

    const nucleotideBaseLabel = svgElementFactory.createTextElement(
      getNucleobaseLabelForCircle(bases[STRAND.SENSE], i), xOfNumbers, Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.SENSE].BASE_LABEL, FONT_SIZES.NUCLEOBASE, determineTextColorBasedOnBackground(bases[STRAND.SENSE][i])
    );

    const linkageStar = phosphorothioateLinkages[STRAND.SENSE][i] ? svgElementFactory.createStarElement(
      baseCircleXCoord + CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS, Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.SENSE].BASE_LABEL + CIRCLE_DIMENSIONS.LINKAGE_STAR_RADIUS, SVG_ELEMENT_COLORS.LINKAGE_STAR
    ) : null;

    const svgElements = [
      nucleotideTextLabel,
      nucleotideCircle,
      nucleotideBaseLabel,
      linkageStar,
    ].filter((element) => element !== null) as SVGElement[];

    image.append(...svgElements);
  }

  image.append(
    phosphorothioateLinkages[STRAND.SENSE][bases[STRAND.SENSE].length] ?
      svgElementFactory.createStarElement(computeNucleobaseCircleXPosition(bases[STRAND.SENSE].length, rightOverhangs[STRAND.SENSE]) +
      CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS, Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.SENSE].BASE_LABEL + CIRCLE_DIMENSIONS.LINKAGE_STAR_RADIUS, SVG_ELEMENT_COLORS.LINKAGE_STAR) : '',
  );

  const numberOfAsNucleotides = bases[STRAND.ANTISENSE].filter((value) => !isOverhangNucleotide(value)).length;
  if (isAntisenseStrandActive) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = bases[STRAND.ANTISENSE].length - 1; i > -1; i--) {
      if (!isOverhangNucleotide(bases[STRAND.ANTISENSE][i]))
        nucleotideCounter--;
      const xOfNumbers = computeNucleobaseCircleXPosition(i, rightOverhangs[STRAND.ANTISENSE]) +
        computeNumericLabelOffset(bases[STRAND.ANTISENSE], i, nucleotideCounter + 1);
      const n = (!isOverhangNucleotide(bases[STRAND.ANTISENSE][i]) && modificationsWithNumericLabels.includes(bases[STRAND.ANTISENSE][i])) ?
        String(nucleotideCounter + 1) : '';
      image.append(
        svgElementFactory.createTextElement(n, xOfNumbers, Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].NUMERIC_LABEL, FONT_SIZES.COMMENT, SVG_ELEMENT_COLORS.TEXT),
        svgElementFactory.createCircleElement(computeNucleobaseCircleXPosition(i, rightOverhangs[STRAND.ANTISENSE]), Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].NUCLEOBASE_CIRCLE, CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS, retrieveNucleobaseColor(bases[STRAND.ANTISENSE][i])),
        svgElementFactory.createTextElement(getNucleobaseLabelForCircle(bases[STRAND.ANTISENSE], i),
          computeNucleobaseCircleXPosition(i, rightOverhangs[STRAND.ANTISENSE]) + computeNumericLabelOffset(bases[STRAND.ANTISENSE], i, nucleotideCounter + 1),
          Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].BASE_LABEL, FONT_SIZES.NUCLEOBASE, determineTextColorBasedOnBackground(bases[STRAND.ANTISENSE][i])),
        phosphorothioateLinkages[STRAND.ANTISENSE][i] ? svgElementFactory.createStarElement(computeNucleobaseCircleXPosition(i, rightOverhangs[STRAND.ANTISENSE]) +
          CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS, Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].BASE_LABEL + CIRCLE_DIMENSIONS.LINKAGE_STAR_RADIUS, SVG_ELEMENT_COLORS.LINKAGE_STAR) : '',
      );
    }
    image.append(
      phosphorothioateLinkages[STRAND.ANTISENSE][bases[STRAND.ANTISENSE].length] ?
        svgElementFactory.createStarElement(computeNucleobaseCircleXPosition(bases[STRAND.ANTISENSE].length, rightOverhangs[STRAND.ANTISENSE]) + CIRCLE_DIMENSIONS.NUCLEOBASE_RADIUS, Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].BASE_LABEL + CIRCLE_DIMENSIONS.LINKAGE_STAR_RADIUS,
          SVG_ELEMENT_COLORS.LINKAGE_STAR) : '',
    );
  }

  const title = `${patternName} for ${numberOfSsNucleotides}${(isAntisenseStrandActive ? `/${numberOfAsNucleotides}` : '')}mer`;
  image.append(svgElementFactory.createTextElement(title, TITLE_POSITION.X, TITLE_POSITION.Y, FONT_SIZES.NUCLEOBASE, SVG_ELEMENT_COLORS.TITLE_TEXT));
  for (let i = 0; i < distinctNucleobaseTypes.length; i++) {
    image.append(
      svgElementFactory.createCircleElement(computeLegendCircleXPosition(i), Y_POSITIONS_FOR_SVG_ELEMENTS.LEGEND_CIRCLES(isAntisenseStrandActive), CIRCLE_DIMENSIONS.LEGEND_RADIUS, retrieveNucleobaseColor(distinctNucleobaseTypes[i])),
      svgElementFactory.createTextElement(distinctNucleobaseTypes[i], computeLegendCircleXPosition(i) + CIRCLE_DIMENSIONS.LEGEND_RADIUS + 4, Y_POSITIONS_FOR_SVG_ELEMENTS.LEGEND_TEXT(isAntisenseStrandActive), FONT_SIZES.COMMENT,
        SVG_ELEMENT_COLORS.TEXT),
    );
  }
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

function getNucleobaseLabelForCircle(bases: string[], index: number): string {
  const base = bases[index];
  const isValidBase = !isOverhangNucleotide(base) && NUCLEOTIDES.includes(base);

  return isValidBase ? base : '';
}

function determineTextColorBasedOnBackground(base: string): string {
  const baseColor = styleMap[base]?.color || '';

  const rgbValues = baseColor.match(/\d+/g)?.map(Number);
  if (!rgbValues || rgbValues.length < 3) {
    return TEXT_COLOR.LIGHT;
  }

  const [r, g, b] = rgbValues;
  const luminance = r * LUMINANCE.RED_COEFFICIENT + g * LUMINANCE.GREEN_COEFFICIENT + b * LUMINANCE.BLUE_COEFFICIENT;
  return luminance > LUMINANCE.THRESHOLD ? TEXT_COLOR.DARK : TEXT_COLOR.LIGHT;
}

function retrieveNucleobaseColor(nucleobase: string): string {
  return styleMap[nucleobase].color;
}

function computeMaxWidthForStrandEnd(strandEnd: STRAND_END): number {
  return Math.max(
    ...STRANDS.map(strand =>
      computeTextWidthInPixels(
        STRAND_END_LABEL_TEXT[strandEnd][strand],
        FONT_SIZES.NUCLEOBASE,
        DEFAULT_FONT_FAMILY
      )
    )
  );
}
