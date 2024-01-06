import {NUCLEOTIDES} from '../../common/model/const';
import {axolabsStyleMap} from '../../common/data-loading-utils/json-loader';
import {SVGElementFactory} from '../view/svg-element-factory';

// dimensions
const NUCLEOBASE_CIRCLE_RADIUS = 15;
const NUCLEOBASE_CIRCLE_DIAMETER = 2 * NUCLEOBASE_CIRCLE_RADIUS;
const LEGEND_CIRCLE_RADIUS = 6;
const LINKAGE_STAR_RADIUS = 5;
const NUCLEOBASE_FONT_SIZE = 17;
const COMMENT_FONT_SIZE = 14;

const COLORS = {
  LINKAGE_STAR: 'red',
  TEXT: 'var(--grey-6)',
  TITLE_TEXT: 'black',
  MODIFICATION_TEXT: 'red'
} as const;

const STRAND_LABELS = {
  LEFT: {
    SENSE_STRAND: 'SS: 5\'',
    ANTISENSE_STRAND: 'AS: 3\'',
  },
  RIGHT: {
    SENSE_STRAND: '3\'',
    ANTISENSE_STRAND: '5\''
  }
} as const;

const NUMERIC_LABEL_SHIFT = {
  ONE_DIGIT: -5,
  TWO_DIGIT: -10,
} as const;

const DEFAULT_FONT_FAMILY = 'Arial';

const STRAND_LABEL_WIDTH = Math.max(
  measureTextWidth(STRAND_LABELS.LEFT.SENSE_STRAND, NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  measureTextWidth(STRAND_LABELS.LEFT.ANTISENSE_STRAND, NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
);

const TERMINUS_LABEL_WIDTH = Math.max(
  measureTextWidth(STRAND_LABELS.RIGHT.SENSE_STRAND, NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
  measureTextWidth(STRAND_LABELS.RIGHT.ANTISENSE_STRAND, NUCLEOBASE_FONT_SIZE, DEFAULT_FONT_FAMILY),
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

const SVG_X_COORDS = {
  LEFT_LABELS: 0,
};

const X_OF_LEFT_MODIFICATIONS = SVG_X_COORDS.LEFT_LABELS + STRAND_LABEL_WIDTH - 5;

const SVG_Y_COORDS = {
  COMMENT: (asExists: boolean) => (asExists) ? 11 * NUCLEOBASE_CIRCLE_RADIUS : 8.5 * NUCLEOBASE_CIRCLE_RADIUS,
  LEGEND_CIRCLES: (asExists: boolean) => (asExists) ? 9.5 * NUCLEOBASE_CIRCLE_RADIUS : 6 * NUCLEOBASE_CIRCLE_RADIUS,
  LEGEND_TEXT: (asExists: boolean) => (asExists) ? 10 * NUCLEOBASE_CIRCLE_RADIUS - 3 : STRAND_DETAILS.ANTISENSE.CIRCLE_Y - 3,
  SVG_TOTAL_HEIGHT: (asExists: boolean) => (asExists) ? 11 * NUCLEOBASE_CIRCLE_RADIUS : 9 * NUCLEOBASE_CIRCLE_RADIUS,
};


export function renderNucleotidePattern(
  patternName: string, antisenseStrandExists: boolean, ssBases: string[],
  asBases: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[],
  ss3PrimeModification: string, ss5PrimeModification: string,
  as3PrimeModification: string, as5PrimeModification: string, comment: string,
  modifiableBases: string[],
): Element {
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

  const distinctBaseTypes = antisenseStrandExists ?
    [...new Set(ssBases.concat(asBases))] :
    [...new Set(ssBases)];

  const ptoLinkageExists = antisenseStrandExists ?
    ssPtoStatuses.concat(asPtoStatuses).includes(true) :
    ssPtoStatuses.includes(true);

  const legendStartIndex = ptoLinkageExists ? 1 : 0;

  const xOfSsRightModifications = ssRightOverhangs * NUCLEOBASE_CIRCLE_DIAMETER + calculateBaseCircleXCoord(-0.5, 0);
  const xOfAsRightModifications = asRightOverhangs * NUCLEOBASE_CIRCLE_DIAMETER + calculateBaseCircleXCoord(-0.5, 0);

  const xOfRightTexts = Math.max(xOfSsRightModifications, xOfAsRightModifications) + widthOfLeftModification +
    NUCLEOBASE_CIRCLE_DIAMETER * widthOfRightOverhangs;

  const width = STRAND_LABEL_WIDTH + widthOfLeftModification + widthOfBases + widthOfRightModification + TERMINUS_LABEL_WIDTH + NUCLEOBASE_CIRCLE_DIAMETER;

  const svgElementFactory = new SVGElementFactory();

  const image = svgElementFactory.createCanvas(width, SVG_Y_COORDS.SVG_TOTAL_HEIGHT(antisenseStrandExists));

  const leftSenseStrandLabel = svgElementFactory.createTextElement(STRAND_LABELS.LEFT.SENSE_STRAND, SVG_X_COORDS.LEFT_LABELS, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT);

  const rightSenseStrandLabel = svgElementFactory.createTextElement(STRAND_LABELS.RIGHT.SENSE_STRAND, xOfRightTexts, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT);

  const leftAntisenseStrandLabel = antisenseStrandExists ? svgElementFactory.createTextElement(STRAND_LABELS.LEFT.ANTISENSE_STRAND, SVG_X_COORDS.LEFT_LABELS, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT) : '';

  const rightAntisenseStrandLabel = antisenseStrandExists ? svgElementFactory.createTextElement(STRAND_LABELS.RIGHT.ANTISENSE_STRAND, xOfRightTexts, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.TEXT) : '';

  const ss5PrimeModificationLabel = svgElementFactory.createTextElement(ss5PrimeModification, X_OF_LEFT_MODIFICATIONS, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT);

  const as3PrimeModificationLabel = antisenseStrandExists ? svgElementFactory.createTextElement(as3PrimeModification, X_OF_LEFT_MODIFICATIONS, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT) : '';

  const ss3PrimeModificationLabel = svgElementFactory.createTextElement(ss3PrimeModification, xOfSsRightModifications, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT);

  const as5PrimeModificationLabel = antisenseStrandExists ? svgElementFactory.createTextElement(as5PrimeModification, xOfAsRightModifications, STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, COLORS.MODIFICATION_TEXT) : '';

  const commentLabel = svgElementFactory.createTextElement(comment, SVG_X_COORDS.LEFT_LABELS, SVG_Y_COORDS.COMMENT(antisenseStrandExists), COMMENT_FONT_SIZE, COLORS.TEXT);

  const starElement = ptoLinkageExists ? svgElementFactory.createStarElement(NUCLEOBASE_CIRCLE_RADIUS, SVG_Y_COORDS.LEGEND_CIRCLES(antisenseStrandExists), COLORS.LINKAGE_STAR) : '';

  const psLinkageLabel = ptoLinkageExists ? svgElementFactory.createTextElement('ps linkage', 2 * NUCLEOBASE_CIRCLE_RADIUS - 8, SVG_Y_COORDS.LEGEND_TEXT(antisenseStrandExists), COMMENT_FONT_SIZE, COLORS.TEXT) : '';

  const elementsToAppend = [
    leftSenseStrandLabel,
    leftAntisenseStrandLabel,
    rightSenseStrandLabel,
    rightAntisenseStrandLabel,
    ss5PrimeModificationLabel,
    as3PrimeModificationLabel,
    ss3PrimeModificationLabel,
    as5PrimeModificationLabel,
    commentLabel,
    starElement,
    psLinkageLabel
  ]

  image.append(...elementsToAppend);

  const numberOfSsNucleotides = ssBases.filter((value) => !isOverhangNucleotide(value)).length;
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = ssBases.length - 1; i > -1; i--) {
    const xOfNumbers = calculateBaseCircleXCoord(i, ssRightOverhangs) +
      calculateNumberLabelShift(ssBases, ssBases.length - i, numberOfSsNucleotides - nucleotideCounter);
    if (!isOverhangNucleotide(ssBases[i]))
      nucleotideCounter--;
    const nucleotideLabel = (!isOverhangNucleotide(ssBases[i]) && modifiableBases.includes(ssBases[i])) ?
      String(numberOfSsNucleotides - nucleotideCounter) : '';
    image.append(
      svgElementFactory.createTextElement(nucleotideLabel, xOfNumbers, STRAND_DETAILS.SENSE.INDEX_Y, COMMENT_FONT_SIZE, COLORS.TEXT),
      svgElementFactory.createCircleElement(calculateBaseCircleXCoord(i, ssRightOverhangs), STRAND_DETAILS.SENSE.CIRCLE_Y, NUCLEOBASE_CIRCLE_RADIUS, getBaseColor(ssBases[i])),
      svgElementFactory.createTextElement(getBaseTextForCircle(ssBases, i), xOfNumbers, STRAND_DETAILS.SENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE,
        determineOptimalTextColor(ssBases[i])),
      ssPtoStatuses[i] ?
        svgElementFactory.createStarElement(calculateBaseCircleXCoord(i, ssRightOverhangs) + NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.SENSE.LABEL_Y + LINKAGE_STAR_RADIUS, COLORS.LINKAGE_STAR) :
        '',
    );
  }
  image.append(
    ssPtoStatuses[ssBases.length] ?
      svgElementFactory.createStarElement(calculateBaseCircleXCoord(ssBases.length, ssRightOverhangs) +
      NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.SENSE.LABEL_Y + LINKAGE_STAR_RADIUS, COLORS.LINKAGE_STAR) : '',
  );

  const numberOfAsNucleotides = asBases.filter((value) => !isOverhangNucleotide(value)).length;
  if (antisenseStrandExists) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = asBases.length - 1; i > -1; i--) {
      if (!isOverhangNucleotide(asBases[i]))
        nucleotideCounter--;
      const xOfNumbers = calculateBaseCircleXCoord(i, asRightOverhangs) +
        calculateNumberLabelShift(asBases, i, nucleotideCounter + 1);
      const n = (!isOverhangNucleotide(asBases[i]) && modifiableBases.includes(asBases[i])) ?
        String(nucleotideCounter + 1) : '';
      image.append(
        svgElementFactory.createTextElement(n, xOfNumbers, STRAND_DETAILS.ANTISENSE.INDEX_Y, COMMENT_FONT_SIZE, COLORS.TEXT),
        svgElementFactory.createCircleElement(calculateBaseCircleXCoord(i, asRightOverhangs), STRAND_DETAILS.ANTISENSE.CIRCLE_Y, NUCLEOBASE_CIRCLE_RADIUS, getBaseColor(asBases[i])),
        svgElementFactory.createTextElement(getBaseTextForCircle(asBases, i),
          calculateBaseCircleXCoord(i, asRightOverhangs) + calculateNumberLabelShift(asBases, i, nucleotideCounter + 1),
          STRAND_DETAILS.ANTISENSE.LABEL_Y, NUCLEOBASE_FONT_SIZE, determineOptimalTextColor(asBases[i])),
        asPtoStatuses[i] ? svgElementFactory.createStarElement(calculateBaseCircleXCoord(i, asRightOverhangs) +
          NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.ANTISENSE.LABEL_Y + LINKAGE_STAR_RADIUS, COLORS.LINKAGE_STAR) : '',
      );
    }
    image.append(
      asPtoStatuses[asBases.length] ?
        svgElementFactory.createStarElement(calculateBaseCircleXCoord(asBases.length, asRightOverhangs) + NUCLEOBASE_CIRCLE_RADIUS, STRAND_DETAILS.ANTISENSE.LABEL_Y + LINKAGE_STAR_RADIUS,
          COLORS.LINKAGE_STAR) : '',
    );
  }

  const title = `${patternName} for ${numberOfSsNucleotides}${(antisenseStrandExists ? `/${numberOfAsNucleotides}` : '')}mer`;
  image.append(svgElementFactory.createTextElement(title, TITLE_POSITION.X, TITLE_POSITION.Y, NUCLEOBASE_FONT_SIZE, COLORS.TITLE_TEXT));
  for (let i = 0; i < distinctBaseTypes.length; i++) {
    image.append(
      svgElementFactory.createCircleElement(calculateLegendXCoord(i), SVG_Y_COORDS.LEGEND_CIRCLES(antisenseStrandExists), LEGEND_CIRCLE_RADIUS, getBaseColor(distinctBaseTypes[i])),
      svgElementFactory.createTextElement(distinctBaseTypes[i], calculateLegendXCoord(i) + LEGEND_CIRCLE_RADIUS + 4, SVG_Y_COORDS.LEGEND_TEXT(antisenseStrandExists), COMMENT_FONT_SIZE,
        COLORS.TEXT),
    );
  }
  return image;
}

function isOverhangNucleotide(modification: string): boolean {
  const overhangSuffix = '(o)';
  return modification.endsWith(overhangSuffix);
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

function determineOptimalTextColor(base: string): string {
  const RED_COEFFICIENT = 0.299;
  const GREEN_COEFFICIENT = 0.587;
  const BLUE_COEFFICIENT = 0.114;
  const LUMINANCE_THRESHOLD = 186;
  const DARK_COLOR = '#333333';
  const LIGHT_COLOR = '#ffffff';

  const styleMap = axolabsStyleMap;
  const baseColor = styleMap[base]?.color || '';
  
  const rgbValues = baseColor.match(/\d+/g)?.map(Number);
  if (!rgbValues || rgbValues.length < 3) {
    return LIGHT_COLOR;
  }

  const [r, g, b] = rgbValues;
  const luminance = r * RED_COEFFICIENT + g * GREEN_COEFFICIENT + b * BLUE_COEFFICIENT;
  return luminance > LUMINANCE_THRESHOLD ? DARK_COLOR : LIGHT_COLOR;
}

function getBaseColor(base: string): string {
  const styleMap = axolabsStyleMap;
  return styleMap[base].color;
}
