import {NUCLEOTIDES} from '../../common/model/const';
import {axolabsStyleMap as styleMap} from '../../common/data-loading-utils/json-loader';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../model/helpers';
import {
  SENSE_STRAND, ANTISENSE_STRAND, STRANDS, THREE_PRIME, FIVE_PRIME, TERMINI
} from '../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../model/types';

const LEFT_END = 'LEFT';
const RIGHT_END = 'RIGHT';
const ENDS = [LEFT_END, RIGHT_END] as const;

// luminance
const RED_COEFFICIENT = 0.299;
const GREEN_COEFFICIENT = 0.587;
const BLUE_COEFFICIENT = 0.114;
const LUMINANCE_THRESHOLD = 186;
const DARK_COLOR = '#333333';
const LIGHT_COLOR = '#ffffff';

// dimensions
const NUCLEOBASE_CIRCLE_RADIUS = 15;
const NUCLEOBASE_CIRCLE_DIAMETER = 2 * NUCLEOBASE_CIRCLE_RADIUS;
const LEGEND_CIRCLE_RADIUS = 6;
const LINKAGE_STAR_RADIUS = 5;
const NUCLEOBASE_FONT_SIZE = 17;
const COMMENT_FONT_SIZE = 14;

// style
const COLORS = {
  LINKAGE_STAR: 'red',
  TEXT: 'var(--grey-6)',
  TITLE_TEXT: 'black',
  MODIFICATION_TEXT: 'red'
} as const;

const STRAND_LABELS = {
  [LEFT_END]: {
    [SENSE_STRAND]: `${SENSE_STRAND}: ${FIVE_PRIME}`,
    [ANTISENSE_STRAND]: `${ANTISENSE_STRAND}: ${THREE_PRIME}`,
  },
  [RIGHT_END]: {
    [SENSE_STRAND]: `${THREE_PRIME}`,
    [ANTISENSE_STRAND]: `${FIVE_PRIME}`,
  }
} as const;

const NUMERIC_LABEL_SHIFT = {
  ONE_DIGIT: -5,
  TWO_DIGIT: -10,
} as const;

const DEFAULT_FONT_FAMILY = 'Arial';

const STRAND_LABEL_WIDTH = Math.max(
  measureTextWidth(STRAND_LABELS[LEFT_END][SENSE_STRAND], NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  measureTextWidth(STRAND_LABELS[LEFT_END][ANTISENSE_STRAND], NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
);

const TERMINUS_LABEL_WIDTH = Math.max(
  measureTextWidth(STRAND_LABELS[RIGHT_END][SENSE_STRAND], NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  measureTextWidth(STRAND_LABELS[RIGHT_END][ANTISENSE_STRAND], NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
);

const TITLE_POSITION = {
  X: NUCLEOBASE_CIRCLE_RADIUS,
  Y: NUCLEOBASE_CIRCLE_RADIUS,
};

const STRAND_DETAILS = {
  SENSE: {
    INDEX_Y: 2 * NUCLEOBASE_CIRCLE_RADIUS,
    CIRCLE_Y: 3.5 * NUCLEOBASE_CIRCLE_RADIUS,
    LABEL_Y: 4 * NUCLEOBASE_CIRCLE_RADIUS,
  },
  ANTISENSE: {
    CIRCLE_Y: 6.5 * NUCLEOBASE_CIRCLE_RADIUS,
    LABEL_Y: 7 * NUCLEOBASE_CIRCLE_RADIUS,
    INDEX_Y: 8.5 * NUCLEOBASE_CIRCLE_RADIUS,
  }
};

const LEFT_LABELS_X_COORDINATE = 0;

const X_OF_LEFT_MODIFICATIONS = LEFT_LABELS_X_COORDINATE + STRAND_LABEL_WIDTH - 5;

const SVG_Y_COORDS = {
  COMMENT: (asExists: boolean) => (asExists) ? 11 * NUCLEOBASE_CIRCLE_RADIUS : 8.5 * NUCLEOBASE_CIRCLE_RADIUS,
  LEGEND_CIRCLES: (asExists: boolean) => (asExists) ? 9.5 * NUCLEOBASE_CIRCLE_RADIUS : 6 * NUCLEOBASE_CIRCLE_RADIUS,
  LEGEND_TEXT: (asExists: boolean) => (asExists) ? 10 * NUCLEOBASE_CIRCLE_RADIUS - 3 : STRAND_DETAILS.ANTISENSE.CIRCLE_Y - 3,
  SVG_TOTAL_HEIGHT: (asExists: boolean) => (asExists) ? 11 * NUCLEOBASE_CIRCLE_RADIUS : 9 * NUCLEOBASE_CIRCLE_RADIUS,
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

  function calculateLegendXCoord(index: number): number {
    const totalPositions = distinctBaseTypes.length + legendStartIndex;
    const spacingUnit = width / totalPositions;
    const position = (index + legendStartIndex) * spacingUnit;
    const adjustedPosition = position + LEGEND_CIRCLE_RADIUS;
    return Math.round(adjustedPosition);
  }

  function calculateBaseCircleXCoord(index: number, rightOverhangs: number): number {
    const rightModificationOffset = maxWidthOfTerminalLabelsByEnd[RIGHT_END];
    const positionalIndex = maxStrandLength - index + rightOverhangs + 1;
    const xCoordinate = positionalIndex * NUCLEOBASE_CIRCLE_DIAMETER;
    const finalPosition = rightModificationOffset + xCoordinate;
    return finalPosition;
  }

  function calculateNumberLabelShift(bases: string[], generalIndex: number, nucleotideIndex: number): number {
    const isSingleDigitOrValidNucleotide = isSingleDigitNumber(nucleotideIndex) || NUCLEOTIDES.includes(bases[generalIndex]);
    const shiftAmount = isSingleDigitOrValidNucleotide 
      ? NUMERIC_LABEL_SHIFT.ONE_DIGIT 
      : NUMERIC_LABEL_SHIFT.TWO_DIGIT;
    return shiftAmount;
  }

  bases[SENSE_STRAND] = bases[SENSE_STRAND].reverse();
  phosphorothioateLinkages[SENSE_STRAND] = phosphorothioateLinkages[SENSE_STRAND].reverse();

  const rightOverhangs = STRANDS.reduce((acc, strand) => {
    acc[strand] = countRightEdgeOverhangNucleotides(bases[strand]);
    return acc;
  }, {} as Record<typeof STRANDS[number], number>);

  const maxStrandLength = Math.max(...STRANDS.map(strand => bases[strand].length - rightOverhangs[strand]));

  const widthOfRightOverhangs = Math.max(rightOverhangs[SENSE_STRAND], rightOverhangs[ANTISENSE_STRAND]);
  const widthOfBases = NUCLEOBASE_CIRCLE_DIAMETER * (maxStrandLength + widthOfRightOverhangs);

  // const widthOfTerminalModification = {} as Record<typeof ENDS[number], number>;
  // widthOfTerminalModification[LEFT_END] = Math.max(
  //   measureTextWidth(terminalModifications[SENSE_STRAND][THREE_PRIME], NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  //   measureTextWidth(terminalModifications[ANTISENSE_STRAND][FIVE_PRIME], NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  // );
  // widthOfTerminalModification[RIGHT_END] = Math.max(
  //   measureTextWidth(terminalModifications[ANTISENSE_STRAND][THREE_PRIME], NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  //   measureTextWidth(terminalModifications[SENSE_STRAND][FIVE_PRIME], NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  // );


  const strandEndToTerminusMapping = {
    [SENSE_STRAND]: {
      [LEFT_END]: THREE_PRIME,
      [RIGHT_END]: FIVE_PRIME
    },
    [ANTISENSE_STRAND]: {
      [LEFT_END]: FIVE_PRIME,
      [RIGHT_END]: THREE_PRIME
    }
  } as const;

  function calculateMaxWidthOfTerminalLabels(end: typeof ENDS[number]) {
    return Math.max(
      ...STRANDS.map(strand => 
        measureTextWidth(
          terminalModifications[strand][strandEndToTerminusMapping[strand][end]],
          NUCLEOBASE_FONT_SIZE, 
          DEFAULT_FONT_FAMILY
        )
      )
    );
  };

  const maxWidthOfTerminalLabelsByEnd = ENDS.reduce((acc, end) => {
    acc[end] = calculateMaxWidthOfTerminalLabels(end);
    return acc;
  }, {} as Record<typeof ENDS[number], number>);

  const distinctBaseTypes = [...new Set(
    bases[SENSE_STRAND].concat(
      isAntisenseStrandActive ? bases[ANTISENSE_STRAND] : []
    )
  )];

  const ptoLinkageExists = isAntisenseStrandActive
    ? [
        ...phosphorothioateLinkages[SENSE_STRAND],
        ...phosphorothioateLinkages[ANTISENSE_STRAND]
      ].some(linkage => linkage)
    : phosphorothioateLinkages[SENSE_STRAND].some(linkage => linkage);

  const legendStartIndex = ptoLinkageExists ? 1 : 0;

  const xOfSsRightModifications = rightOverhangs[SENSE_STRAND] * NUCLEOBASE_CIRCLE_DIAMETER + calculateBaseCircleXCoord(-0.5, 0);
  const xOfAsRightModifications = rightOverhangs[ANTISENSE_STRAND] * NUCLEOBASE_CIRCLE_DIAMETER + calculateBaseCircleXCoord(-0.5, 0);

  const rightLabelsXCoordinate = Math.max(xOfSsRightModifications, xOfAsRightModifications) + maxWidthOfTerminalLabelsByEnd[LEFT_END] +
    NUCLEOBASE_CIRCLE_DIAMETER * widthOfRightOverhangs;

  const width = STRAND_LABEL_WIDTH + maxWidthOfTerminalLabelsByEnd[LEFT_END] + widthOfBases + maxWidthOfTerminalLabelsByEnd[RIGHT_END] + TERMINUS_LABEL_WIDTH + NUCLEOBASE_CIRCLE_DIAMETER;

  const svgFactory = new SVGElementFactory();

  const image = svgFactory.createCanvas(width, SVG_Y_COORDS.SVG_TOTAL_HEIGHT(isAntisenseStrandActive));

  const strandLabel = {} as Record<typeof STRANDS[number], Record<typeof ENDS[number], SVGElement | null>>;

  STRANDS.forEach((strand) => {
    strandLabel[strand] = {} as Record<typeof ENDS[number], SVGElement | null>;
    ENDS.forEach((end) => {
      strandLabel[strand][end] = null;
    });
  });

  strandLabel[SENSE_STRAND][LEFT_END]
    = svgFactory.createTextElement(STRAND_LABELS[LEFT_END][SENSE_STRAND], LEFT_LABELS_X_COORDINATE, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT);

  strandLabel[SENSE_STRAND][RIGHT_END]
    = svgFactory.createTextElement(STRAND_LABELS[RIGHT_END][SENSE_STRAND], rightLabelsXCoordinate, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT);

  strandLabel[ANTISENSE_STRAND][LEFT_END]
    = isAntisenseStrandActive ? svgFactory.createTextElement(STRAND_LABELS[LEFT_END][ANTISENSE_STRAND], LEFT_LABELS_X_COORDINATE, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT) : null;

  strandLabel[ANTISENSE_STRAND][RIGHT_END]
    = isAntisenseStrandActive ? svgFactory.createTextElement(STRAND_LABELS[RIGHT_END][ANTISENSE_STRAND], rightLabelsXCoordinate, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT) : null;

  const modificationLabel = {} as Record<typeof STRANDS[number], Record<typeof TERMINI[number], SVGElement | null>>;

  STRANDS.forEach((strand) => {
    modificationLabel[strand] = {} as Record<typeof TERMINI[number], SVGElement | null>;
    TERMINI.forEach((terminus) => {
      modificationLabel[strand][terminus] = null;
    });
  });


  modificationLabel[SENSE_STRAND][FIVE_PRIME]
    = svgFactory.createTextElement(terminalModifications[SENSE_STRAND][FIVE_PRIME], X_OF_LEFT_MODIFICATIONS, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT);

  modificationLabel[ANTISENSE_STRAND][THREE_PRIME]
    = isAntisenseStrandActive ? svgFactory.createTextElement(terminalModifications[ANTISENSE_STRAND][THREE_PRIME], X_OF_LEFT_MODIFICATIONS, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT) : null;

  modificationLabel[SENSE_STRAND][THREE_PRIME]
    = svgFactory.createTextElement(terminalModifications[SENSE_STRAND][THREE_PRIME], xOfSsRightModifications, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT);

  modificationLabel[ANTISENSE_STRAND][FIVE_PRIME]
    = isAntisenseStrandActive ? svgFactory.createTextElement(terminalModifications[ANTISENSE_STRAND][FIVE_PRIME], xOfAsRightModifications, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT) : null;

  const commentLabel = svgFactory.createTextElement(comment, LEFT_LABELS_X_COORDINATE, SVG_Y_COORDS.COMMENT(isAntisenseStrandActive), COMMENT_FONT_SIZE, COLORS.TEXT);

  const starElement = ptoLinkageExists ? svgFactory.createStarElement(NUCLEOBASE_CIRCLE_RADIUS, SVG_Y_COORDS.LEGEND_CIRCLES(isAntisenseStrandActive), COLORS.LINKAGE_STAR) : null;

  const psLinkageLabel = ptoLinkageExists ? svgFactory.createTextElement('ps linkage', 2 * NUCLEOBASE_CIRCLE_RADIUS - 8, SVG_Y_COORDS.LEGEND_TEXT(isAntisenseStrandActive), COMMENT_FONT_SIZE, COLORS.TEXT) : null;

  const elementsToAppend = [
    strandLabel[SENSE_STRAND][LEFT_END],
    strandLabel[SENSE_STRAND][RIGHT_END],
    strandLabel[ANTISENSE_STRAND][LEFT_END],
    strandLabel[ANTISENSE_STRAND][RIGHT_END],
    modificationLabel[SENSE_STRAND][FIVE_PRIME],
    modificationLabel[ANTISENSE_STRAND][THREE_PRIME],
    modificationLabel[SENSE_STRAND][THREE_PRIME],
    modificationLabel[ANTISENSE_STRAND][FIVE_PRIME],
    commentLabel,
    starElement,
    psLinkageLabel,
  ].filter((element) => element !== null) as SVGElement[];

  image.append(...elementsToAppend);

  const numberOfSsNucleotides = bases[SENSE_STRAND].filter((value) => !isOverhangNucleotide(value)).length;
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = bases[SENSE_STRAND].length - 1; i > -1; i--) {
    const xOfNumbers = calculateBaseCircleXCoord(i, rightOverhangs[SENSE_STRAND]) +
      calculateNumberLabelShift(bases[SENSE_STRAND], bases[SENSE_STRAND].length - i, numberOfSsNucleotides - nucleotideCounter);
    if (!isOverhangNucleotide(bases[SENSE_STRAND][i]))
      nucleotideCounter--;
    const nucleotideLabel = (!isOverhangNucleotide(bases[SENSE_STRAND][i]) && modificationsWithNumericLabels.includes(bases[SENSE_STRAND][i])) ?
      String(numberOfSsNucleotides - nucleotideCounter) : '';

    const baseCircleXCoord = calculateBaseCircleXCoord(i, rightOverhangs[SENSE_STRAND]);

    const nucleotideTextLabel = svgFactory.createTextElement(
      nucleotideLabel, xOfNumbers, STRAND_DETAILS.SENSE.INDEX_Y, COMMENT_FONT_SIZE, COLORS.TEXT
    );

    const nucleotideCircle = svgFactory.createCircleElement(
      baseCircleXCoord, STRAND_DETAILS.SENSE.CIRCLE_Y, NUCLEOBASE_CIRCLE_RADIUS, getBaseColor(bases[SENSE_STRAND][i])
    );

    const nucleotideBaseLabel = svgFactory.createTextElement(
      getBaseTextForCircle(bases[SENSE_STRAND], i), xOfNumbers, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, getTextColorForBaseBackground(bases[SENSE_STRAND][i])
    );

    const linkageStar = phosphorothioateLinkages[SENSE_STRAND][i] ? svgFactory.createStarElement(
      baseCircleXCoord + NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.SENSE.LABEL_Y + LINKAGE_STAR_RADIUS, COLORS.LINKAGE_STAR
    ) : null;

    const elementsToAppend = [
      nucleotideTextLabel,
      nucleotideCircle,
      nucleotideBaseLabel,
      linkageStar,
    ].filter((element) => element !== null) as SVGElement[];

    image.append(...elementsToAppend);
  }

  image.append(
    phosphorothioateLinkages[SENSE_STRAND][bases[SENSE_STRAND].length] ?
      svgFactory.createStarElement(calculateBaseCircleXCoord(bases[SENSE_STRAND].length, rightOverhangs[SENSE_STRAND]) +
      NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.SENSE.LABEL_Y + LINKAGE_STAR_RADIUS, COLORS.LINKAGE_STAR) : '',
  );

  const numberOfAsNucleotides = bases[ANTISENSE_STRAND].filter((value) => !isOverhangNucleotide(value)).length;
  if (isAntisenseStrandActive) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = bases[ANTISENSE_STRAND].length - 1; i > -1; i--) {
      if (!isOverhangNucleotide(bases[ANTISENSE_STRAND][i]))
        nucleotideCounter--;
      const xOfNumbers = calculateBaseCircleXCoord(i, rightOverhangs[ANTISENSE_STRAND]) +
        calculateNumberLabelShift(bases[ANTISENSE_STRAND], i, nucleotideCounter + 1);
      const n = (!isOverhangNucleotide(bases[ANTISENSE_STRAND][i]) && modificationsWithNumericLabels.includes(bases[ANTISENSE_STRAND][i])) ?
        String(nucleotideCounter + 1) : '';
      image.append(
        svgFactory.createTextElement(n, xOfNumbers, STRAND_DETAILS.ANTISENSE.INDEX_Y, COMMENT_FONT_SIZE, COLORS.TEXT),
        svgFactory.createCircleElement(calculateBaseCircleXCoord(i, rightOverhangs[ANTISENSE_STRAND]), STRAND_DETAILS.ANTISENSE.CIRCLE_Y, NUCLEOBASE_CIRCLE_RADIUS, getBaseColor(bases[ANTISENSE_STRAND][i])),
        svgFactory.createTextElement(getBaseTextForCircle(bases[ANTISENSE_STRAND], i),
          calculateBaseCircleXCoord(i, rightOverhangs[ANTISENSE_STRAND]) + calculateNumberLabelShift(bases[ANTISENSE_STRAND], i, nucleotideCounter + 1),
          STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, getTextColorForBaseBackground(bases[ANTISENSE_STRAND][i])),
        phosphorothioateLinkages[ANTISENSE_STRAND][i] ? svgFactory.createStarElement(calculateBaseCircleXCoord(i, rightOverhangs[ANTISENSE_STRAND]) +
          NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.ANTISENSE.LABEL_Y + LINKAGE_STAR_RADIUS, COLORS.LINKAGE_STAR) : '',
      );
    }
    image.append(
      phosphorothioateLinkages[ANTISENSE_STRAND][bases[ANTISENSE_STRAND].length] ?
        svgFactory.createStarElement(calculateBaseCircleXCoord(bases[ANTISENSE_STRAND].length, rightOverhangs[ANTISENSE_STRAND]) + NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.ANTISENSE.LABEL_Y + LINKAGE_STAR_RADIUS,
          COLORS.LINKAGE_STAR) : '',
    );
  }

  const title = `${patternName} for ${numberOfSsNucleotides}${(isAntisenseStrandActive ? `/${numberOfAsNucleotides}` : '')}mer`;
  image.append(svgFactory.createTextElement(title, TITLE_POSITION.X, TITLE_POSITION.Y, NUCLEOBASE_FONT_SIZE, COLORS.TITLE_TEXT));
  for (let i = 0; i < distinctBaseTypes.length; i++) {
    image.append(
      svgFactory.createCircleElement(calculateLegendXCoord(i), SVG_Y_COORDS.LEGEND_CIRCLES(isAntisenseStrandActive), LEGEND_CIRCLE_RADIUS, getBaseColor(distinctBaseTypes[i])),
      svgFactory.createTextElement(distinctBaseTypes[i], calculateLegendXCoord(i) + LEGEND_CIRCLE_RADIUS + 4, SVG_Y_COORDS.LEGEND_TEXT(isAntisenseStrandActive), COMMENT_FONT_SIZE,
        COLORS.TEXT),
    );
  }
  return image;
}

function isSingleDigitNumber(n: number): boolean {
  return n >= 0 && n < 10;
}

function countRightEdgeOverhangNucleotides(modifications: string[]): number {
  const lastIdx = modifications.length - 1;
  let count = 0;
  while (count <= lastIdx && isOverhangNucleotide(modifications[count])) {
    count++;
  }
  return count === lastIdx + 1 ? 0 : count;
}

function measureTextWidth(text: string, fontSize: number, fontFamily: string): number {
  const canvas = document.createElement('canvas');
  const context = canvas.getContext('2d');
  if (context) {
    context.font = `${fontSize}px ${fontFamily}`;
    const metrics = context.measureText(text);
    return 2 * metrics.width;
  }
  return 0;
}

function getBaseTextForCircle(bases: string[], index: number): string {
  const base = bases[index];
  const isValidBase = !isOverhangNucleotide(base) && NUCLEOTIDES.includes(base);
  
  return isValidBase ? base : '';
}

function getTextColorForBaseBackground(base: string): string {
  const baseColor = styleMap[base]?.color || '';
  
  const rgbValues = baseColor.match(/\d+/g)?.map(Number);
  if (!rgbValues || rgbValues.length < 3) {
    return LIGHT_COLOR;
  }

  const [r, g, b] = rgbValues;
  const luminance = r * RED_COEFFICIENT + g * GREEN_COEFFICIENT + b * BLUE_COEFFICIENT;
  return luminance > LUMINANCE_THRESHOLD ? DARK_COLOR : LIGHT_COLOR;
}

function getBaseColor(nucleobase: string): string {
  return styleMap[nucleobase].color;
}
