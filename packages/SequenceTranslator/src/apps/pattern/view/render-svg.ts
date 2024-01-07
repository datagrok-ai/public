import {NUCLEOTIDES} from '../../common/model/const';
import {axolabsStyleMap as styleMap} from '../../common/data-loading-utils/json-loader';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../model/helpers';

const SENSE_STRAND = 'SS';
const ANTISENSE_STRAND = 'AS';
const STRANDS = [SENSE_STRAND, ANTISENSE_STRAND] as const;
const FIVE_PRIME = '5\'';
const THREE_PRIME = '3\'';
const TERMINI = [FIVE_PRIME, THREE_PRIME] as const;
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

type bases = Record<typeof STRANDS[number], string>;
type ptoStatuses = Record<typeof STRANDS[number], boolean[]>;
type modifications = Record<typeof STRANDS[number], Record<typeof TERMINI[number], string>>;

export function renderNucleotidePattern(options: {
  patternName: string,
  isAsStrandActive: boolean,
  baseInputValues: Record<typeof STRANDS[number], string[]>,
  ptoLinkageValues: Record<typeof STRANDS[number], boolean[]>,
  terminalModificationValues: Record<typeof STRANDS[number], Record<typeof TERMINI[number], string>>,
  comment: string,
  modificationsWithNumericLabels: string[],
}): Element {
  const {
    patternName,
    isAsStrandActive,
    baseInputValues: {SS: ssBases, AS: asBases},
    ptoLinkageValues: {SS: ssPtoStatuses, AS: asPtoStatuses},
    terminalModificationValues: {[SENSE_STRAND]: {[FIVE_PRIME]: ss5PrimeModification, [THREE_PRIME]: ss3PrimeModification}, [ANTISENSE_STRAND]: {[FIVE_PRIME]: as5PrimeModification, [THREE_PRIME]: as3PrimeModification}},
    comment,
    modificationsWithNumericLabels
  } = options;

  function calculateLegendXCoord(index: number): number {
    const totalPositions = distinctBaseTypes.length + legendStartIndex;
    const spacingUnit = width / totalPositions;
    const position = (index + legendStartIndex) * spacingUnit;
    const adjustedPosition = position + LEGEND_CIRCLE_RADIUS;
    return Math.round(adjustedPosition);
  }

  function calculateBaseCircleXCoord(index: number, rightOverhangs: number): number {
    const rightModificationOffset = widthOfRightModification;
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

  ssBases = ssBases.reverse();
  ssPtoStatuses = ssPtoStatuses.reverse();

  const ssRightOverhangs = countRightEdgeOverhangNucleotides(ssBases);
  const asRightOverhangs = countRightEdgeOverhangNucleotides(asBases);

  const maxStrandLength = Math.max(
    ssBases.length - ssRightOverhangs,
    asBases.length - asRightOverhangs,
  );

  const widthOfRightOverhangs = Math.max(ssRightOverhangs, asRightOverhangs);
  const widthOfBases = NUCLEOBASE_CIRCLE_DIAMETER * (maxStrandLength + widthOfRightOverhangs);

  const widthOfLeftModification = Math.max(
    measureTextWidth(ss3PrimeModification, NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
    measureTextWidth(as5PrimeModification, NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  );

  const widthOfRightModification = Math.max(
    measureTextWidth(ss5PrimeModification, NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
    measureTextWidth(as3PrimeModification, NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  );

  const distinctBaseTypes = isAsStrandActive ?
    [...new Set(ssBases.concat(asBases))] :
    [...new Set(ssBases)];

  const ptoLinkageExists = isAsStrandActive ?
    ssPtoStatuses.concat(asPtoStatuses).includes(true) :
    ssPtoStatuses.includes(true);

  const legendStartIndex = ptoLinkageExists ? 1 : 0;

  const xOfSsRightModifications = ssRightOverhangs * NUCLEOBASE_CIRCLE_DIAMETER + calculateBaseCircleXCoord(-0.5, 0);
  const xOfAsRightModifications = asRightOverhangs * NUCLEOBASE_CIRCLE_DIAMETER + calculateBaseCircleXCoord(-0.5, 0);

  const rightLabelsXCoordinate = Math.max(xOfSsRightModifications, xOfAsRightModifications) + widthOfLeftModification +
    NUCLEOBASE_CIRCLE_DIAMETER * widthOfRightOverhangs;

  const width = STRAND_LABEL_WIDTH + widthOfLeftModification + widthOfBases + widthOfRightModification + TERMINUS_LABEL_WIDTH + NUCLEOBASE_CIRCLE_DIAMETER;

  const svgFactory = new SVGElementFactory();

  const image = svgFactory.createCanvas(width, SVG_Y_COORDS.SVG_TOTAL_HEIGHT(isAsStrandActive));

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
    = isAsStrandActive ? svgFactory.createTextElement(STRAND_LABELS[LEFT_END][ANTISENSE_STRAND], LEFT_LABELS_X_COORDINATE, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT) : null;

  strandLabel[ANTISENSE_STRAND][RIGHT_END]
    = isAsStrandActive ? svgFactory.createTextElement(STRAND_LABELS[RIGHT_END][ANTISENSE_STRAND], rightLabelsXCoordinate, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT) : null;

  const modificationLabel = {} as Record<typeof STRANDS[number], Record<typeof TERMINI[number], SVGElement | null>>;

  STRANDS.forEach((strand) => {
    modificationLabel[strand] = {} as Record<typeof TERMINI[number], SVGElement | null>;
    TERMINI.forEach((terminus) => {
      modificationLabel[strand][terminus] = null;
    });
  });


  modificationLabel[SENSE_STRAND][FIVE_PRIME]
    = svgFactory.createTextElement(ss5PrimeModification, X_OF_LEFT_MODIFICATIONS, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT);

  modificationLabel[ANTISENSE_STRAND][THREE_PRIME]
    = isAsStrandActive ? svgFactory.createTextElement(as3PrimeModification, X_OF_LEFT_MODIFICATIONS, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT) : null;

  modificationLabel[SENSE_STRAND][THREE_PRIME]
    = svgFactory.createTextElement(ss3PrimeModification, xOfSsRightModifications, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT);

  modificationLabel[ANTISENSE_STRAND][FIVE_PRIME]
    = isAsStrandActive ? svgFactory.createTextElement(as5PrimeModification, xOfAsRightModifications, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT) : null;

  const commentLabel = svgFactory.createTextElement(comment, LEFT_LABELS_X_COORDINATE, SVG_Y_COORDS.COMMENT(isAsStrandActive), COMMENT_FONT_SIZE, COLORS.TEXT);

  const starElement = ptoLinkageExists ? svgFactory.createStarElement(NUCLEOBASE_CIRCLE_RADIUS, SVG_Y_COORDS.LEGEND_CIRCLES(isAsStrandActive), COLORS.LINKAGE_STAR) : null;

  const psLinkageLabel = ptoLinkageExists ? svgFactory.createTextElement('ps linkage', 2 * NUCLEOBASE_CIRCLE_RADIUS - 8, SVG_Y_COORDS.LEGEND_TEXT(isAsStrandActive), COMMENT_FONT_SIZE, COLORS.TEXT) : null;

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

  const numberOfSsNucleotides = ssBases.filter((value) => !isOverhangNucleotide(value)).length;
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = ssBases.length - 1; i > -1; i--) {
    const xOfNumbers = calculateBaseCircleXCoord(i, ssRightOverhangs) +
      calculateNumberLabelShift(ssBases, ssBases.length - i, numberOfSsNucleotides - nucleotideCounter);
    if (!isOverhangNucleotide(ssBases[i]))
      nucleotideCounter--;
    const nucleotideLabel = (!isOverhangNucleotide(ssBases[i]) && modificationsWithNumericLabels.includes(ssBases[i])) ?
      String(numberOfSsNucleotides - nucleotideCounter) : '';

    const baseCircleXCoord = calculateBaseCircleXCoord(i, ssRightOverhangs);

    const nucleotideTextLabel = svgFactory.createTextElement(
      nucleotideLabel, xOfNumbers, STRAND_DETAILS.SENSE.INDEX_Y, COMMENT_FONT_SIZE, COLORS.TEXT
    );

    const nucleotideCircle = svgFactory.createCircleElement(
      baseCircleXCoord, STRAND_DETAILS.SENSE.CIRCLE_Y, NUCLEOBASE_CIRCLE_RADIUS, getBaseColor(ssBases[i])
    );

    const nucleotideBaseLabel = svgFactory.createTextElement(
      getBaseTextForCircle(ssBases, i), xOfNumbers, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, getTextColorForBaseBackground(ssBases[i])
    );

    const linkageStar = ssPtoStatuses[i] ? svgFactory.createStarElement(
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
    ssPtoStatuses[ssBases.length] ?
      svgFactory.createStarElement(calculateBaseCircleXCoord(ssBases.length, ssRightOverhangs) +
      NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.SENSE.LABEL_Y + LINKAGE_STAR_RADIUS, COLORS.LINKAGE_STAR) : '',
  );

  const numberOfAsNucleotides = asBases.filter((value) => !isOverhangNucleotide(value)).length;
  if (isAsStrandActive) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = asBases.length - 1; i > -1; i--) {
      if (!isOverhangNucleotide(asBases[i]))
        nucleotideCounter--;
      const xOfNumbers = calculateBaseCircleXCoord(i, asRightOverhangs) +
        calculateNumberLabelShift(asBases, i, nucleotideCounter + 1);
      const n = (!isOverhangNucleotide(asBases[i]) && modificationsWithNumericLabels.includes(asBases[i])) ?
        String(nucleotideCounter + 1) : '';
      image.append(
        svgFactory.createTextElement(n, xOfNumbers, STRAND_DETAILS.ANTISENSE.INDEX_Y, COMMENT_FONT_SIZE, COLORS.TEXT),
        svgFactory.createCircleElement(calculateBaseCircleXCoord(i, asRightOverhangs), STRAND_DETAILS.ANTISENSE.CIRCLE_Y, NUCLEOBASE_CIRCLE_RADIUS, getBaseColor(asBases[i])),
        svgFactory.createTextElement(getBaseTextForCircle(asBases, i),
          calculateBaseCircleXCoord(i, asRightOverhangs) + calculateNumberLabelShift(asBases, i, nucleotideCounter + 1),
          STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, getTextColorForBaseBackground(asBases[i])),
        asPtoStatuses[i] ? svgFactory.createStarElement(calculateBaseCircleXCoord(i, asRightOverhangs) +
          NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.ANTISENSE.LABEL_Y + LINKAGE_STAR_RADIUS, COLORS.LINKAGE_STAR) : '',
      );
    }
    image.append(
      asPtoStatuses[asBases.length] ?
        svgFactory.createStarElement(calculateBaseCircleXCoord(asBases.length, asRightOverhangs) + NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.ANTISENSE.LABEL_Y + LINKAGE_STAR_RADIUS,
          COLORS.LINKAGE_STAR) : '',
    );
  }

  const title = `${patternName} for ${numberOfSsNucleotides}${(isAsStrandActive ? `/${numberOfAsNucleotides}` : '')}mer`;
  image.append(svgFactory.createTextElement(title, TITLE_POSITION.X, TITLE_POSITION.Y, NUCLEOBASE_FONT_SIZE, COLORS.TITLE_TEXT));
  for (let i = 0; i < distinctBaseTypes.length; i++) {
    image.append(
      svgFactory.createCircleElement(calculateLegendXCoord(i), SVG_Y_COORDS.LEGEND_CIRCLES(isAsStrandActive), LEGEND_CIRCLE_RADIUS, getBaseColor(distinctBaseTypes[i])),
      svgFactory.createTextElement(distinctBaseTypes[i], calculateLegendXCoord(i) + LEGEND_CIRCLE_RADIUS + 4, SVG_Y_COORDS.LEGEND_TEXT(isAsStrandActive), COMMENT_FONT_SIZE,
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
