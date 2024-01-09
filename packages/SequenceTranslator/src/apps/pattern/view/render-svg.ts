import {NUCLEOTIDES} from '../../common/model/const';
import {axolabsStyleMap as styleMap} from '../../common/data-loading-utils/json-loader';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../model/helpers';
import {
  SENSE_STRAND, ANTISENSE_STRAND, STRANDS, THREE_PRIME, FIVE_PRIME, TERMINI
} from '../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../model/types';

const END_LEFT = 'LEFT';
const END_RIGHT = 'RIGHT';
const ENDS = [END_LEFT, END_RIGHT] as const;

// luminance
const RED_COEFFICIENT = 0.299;
const GREEN_COEFFICIENT = 0.587;
const BLUE_COEFFICIENT = 0.114;
const LUMINANCE_THRESHOLD = 186;
const DARK_COLOR = '#333333';
const LIGHT_COLOR = '#ffffff';

// dimensions
const NUCLEOBASE_CIRCLE_RADIUS_PX = 15;
const NUCLEOBASE_CIRCLE_DIAMETER_PX = 2 * NUCLEOBASE_CIRCLE_RADIUS_PX;
const LEGEND_CIRCLE_RADIUS_PX = 6;
const LINKAGE_STAR_RADIUS_PX = 5;
const NUCLEOBASE_FONT_SIZE_PX = 17;
const COMMENT_FONT_SIZE_PX = 14;

// style
const SVG_ELEMENT_COLORS = {
  LINKAGE_STAR: 'red',
  TEXT: 'var(--grey-6)',
  TITLE_TEXT: 'black',
  MODIFICATION_TEXT: 'red'
} as const;

const STRAND_END_LABEL_TEXT = {
  [END_LEFT]: {
    [SENSE_STRAND]: `${SENSE_STRAND}: ${FIVE_PRIME}`,
    [ANTISENSE_STRAND]: `${ANTISENSE_STRAND}: ${THREE_PRIME}`,
  },
  [END_RIGHT]: {
    [SENSE_STRAND]: `${THREE_PRIME}`,
    [ANTISENSE_STRAND]: `${FIVE_PRIME}`,
  }
} as const;

const NUMERIC_LABEL_POSITION_OFFSET = {
  ONE_DIGIT: -5,
  TWO_DIGIT: -10,
} as const;

const DEFAULT_FONT_FAMILY = 'Arial';

const WIDTH_STRAND_LABEL = Math.max(
  computeTextWidthInPixels(STRAND_END_LABEL_TEXT[END_LEFT][SENSE_STRAND], NUCLEOBASE_FONT_SIZE_PX, DEFAULT_FONT_FAMILY),
  computeTextWidthInPixels(STRAND_END_LABEL_TEXT[END_LEFT][ANTISENSE_STRAND], NUCLEOBASE_FONT_SIZE_PX, DEFAULT_FONT_FAMILY),
);

const WIDTH_TERMINUS_LABEL = Math.max(
  computeTextWidthInPixels(STRAND_END_LABEL_TEXT[END_RIGHT][SENSE_STRAND], NUCLEOBASE_FONT_SIZE_PX, DEFAULT_FONT_FAMILY),
  computeTextWidthInPixels(STRAND_END_LABEL_TEXT[END_RIGHT][ANTISENSE_STRAND], NUCLEOBASE_FONT_SIZE_PX, DEFAULT_FONT_FAMILY),
);

const TITLE_POSITION = {
  X: NUCLEOBASE_CIRCLE_RADIUS_PX,
  Y: NUCLEOBASE_CIRCLE_RADIUS_PX,
};

const POSITIONS_FOR_STRAND_ELEMENTS = {
  [SENSE_STRAND]: {
    Y_POSITION_NUMERIC_LABEL: 2 * NUCLEOBASE_CIRCLE_RADIUS_PX,
    Y_POSITION_NUCLEOBASE_CIRCLE: 3.5 * NUCLEOBASE_CIRCLE_RADIUS_PX,
    Y_POSITION_BASE_LABEL: 4 * NUCLEOBASE_CIRCLE_RADIUS_PX,
  },
  [ANTISENSE_STRAND]: {
    Y_POSITION_NUMERIC_LABEL: 8.5 * NUCLEOBASE_CIRCLE_RADIUS_PX,
    Y_POSITION_NUCLEOBASE_CIRCLE: 6.5 * NUCLEOBASE_CIRCLE_RADIUS_PX,
    Y_POSITION_BASE_LABEL: 7 * NUCLEOBASE_CIRCLE_RADIUS_PX,
  }
};

const SVG_VERTICAL_COORDINATES = {
  COMMENT: (asExists: boolean) => (asExists) ? 11 * NUCLEOBASE_CIRCLE_RADIUS_PX : 8.5 * NUCLEOBASE_CIRCLE_RADIUS_PX,
  LEGEND_CIRCLES: (asExists: boolean) => (asExists) ? 9.5 * NUCLEOBASE_CIRCLE_RADIUS_PX : 6 * NUCLEOBASE_CIRCLE_RADIUS_PX,
  LEGEND_TEXT: (asExists: boolean) => (asExists) ? 10 * NUCLEOBASE_CIRCLE_RADIUS_PX - 3 : POSITIONS_FOR_STRAND_ELEMENTS[ANTISENSE_STRAND].Y_POSITION_NUCLEOBASE_CIRCLE - 3,
  SVG_TOTAL_HEIGHT: (asExists: boolean) => (asExists) ? 11 * NUCLEOBASE_CIRCLE_RADIUS_PX : 9 * NUCLEOBASE_CIRCLE_RADIUS_PX,
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
    const adjustedPosition = position + LEGEND_CIRCLE_RADIUS_PX;
    return Math.round(adjustedPosition);
  }

  function computeNucleobaseCircleXPosition(index: number, rightOverhangs: number): number {
    const rightModificationOffset = maxWidthOfTerminusLabelsByEnd[END_RIGHT];
    const positionalIndex = maxStrandLength - index + rightOverhangs + 1;
    const xPosition = positionalIndex * NUCLEOBASE_CIRCLE_DIAMETER_PX;
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

  bases[SENSE_STRAND] = bases[SENSE_STRAND].reverse();
  phosphorothioateLinkages[SENSE_STRAND] = phosphorothioateLinkages[SENSE_STRAND].reverse();

  const rightOverhangs = STRANDS.reduce((acc, strand) => {
    acc[strand] = countOverhangNucleotidesAtStrandEnd(bases[strand]);
    return acc;
  }, {} as Record<typeof STRANDS[number], number>);

  const maxStrandLength = Math.max(...STRANDS.map(strand => bases[strand].length - rightOverhangs[strand]));

  const widthOfRightOverhangs = Math.max(rightOverhangs[SENSE_STRAND], rightOverhangs[ANTISENSE_STRAND]);
  const widthOfBases = NUCLEOBASE_CIRCLE_DIAMETER_PX * (maxStrandLength + widthOfRightOverhangs);

  const strandEndToTerminusMapping = {
    [SENSE_STRAND]: {
      [END_LEFT]: THREE_PRIME,
      [END_RIGHT]: FIVE_PRIME
    },
    [ANTISENSE_STRAND]: {
      [END_LEFT]: FIVE_PRIME,
      [END_RIGHT]: THREE_PRIME
    }
  } as const;

  function computeMaxWidthOfTerminalLabels(end: typeof ENDS[number]) {
    return Math.max(
      ...STRANDS.map(strand =>
        computeTextWidthInPixels(
          terminalModifications[strand][strandEndToTerminusMapping[strand][end]],
          NUCLEOBASE_FONT_SIZE_PX,
          DEFAULT_FONT_FAMILY
        )
      )
    );
  };

  const maxWidthOfTerminusLabelsByEnd = ENDS.reduce((acc, end) => {
    acc[end] = computeMaxWidthOfTerminalLabels(end);
    return acc;
  }, {} as Record<typeof ENDS[number], number>);

  const distinctNucleobaseTypes = [...new Set(
    bases[SENSE_STRAND].concat(
      isAntisenseStrandActive ? bases[ANTISENSE_STRAND] : []
    )
  )];

  const isPhosphorothioateLinkageActive = isAntisenseStrandActive
    ? [
        ...phosphorothioateLinkages[SENSE_STRAND],
        ...phosphorothioateLinkages[ANTISENSE_STRAND]
      ].some(linkage => linkage)
    : phosphorothioateLinkages[SENSE_STRAND].some(linkage => linkage);

  const legendStartIndex = isPhosphorothioateLinkageActive ? 1 : 0;

  const xPositionOfTerminusModifications = {} as Record<typeof ENDS[number], Record<typeof STRANDS[number], number>>;
  xPositionOfTerminusModifications[END_LEFT] = STRANDS.reduce((acc, strand) => {
    acc[strand] = WIDTH_STRAND_LABEL - 5;
    return acc;
  }, {} as Record<typeof STRANDS[number], number>);
  xPositionOfTerminusModifications[END_RIGHT] = STRANDS.reduce((acc, strand) => {
    acc[strand] = rightOverhangs[strand] * NUCLEOBASE_CIRCLE_DIAMETER_PX + computeNucleobaseCircleXPosition(-0.5, 0);
    return acc;
  }, {} as Record<typeof STRANDS[number], number>);

  const xPositionOfStrandLabels = {
    [END_LEFT]: 0,
    [END_RIGHT]: Math.max(...STRANDS.map(strand => xPositionOfTerminusModifications[END_RIGHT][strand]))
      + maxWidthOfTerminusLabelsByEnd[END_LEFT]
      + NUCLEOBASE_CIRCLE_DIAMETER_PX * widthOfRightOverhangs
  };

  const width = WIDTH_STRAND_LABEL + maxWidthOfTerminusLabelsByEnd[END_LEFT] + widthOfBases + maxWidthOfTerminusLabelsByEnd[END_RIGHT] + WIDTH_TERMINUS_LABEL + NUCLEOBASE_CIRCLE_DIAMETER_PX;

  const svgElementFactory = new SVGElementFactory();

  const image = svgElementFactory.createCanvas(width, SVG_VERTICAL_COORDINATES.SVG_TOTAL_HEIGHT(isAntisenseStrandActive));

  const labelsStrandEnd = {} as Record<typeof STRANDS[number], Record<typeof ENDS[number], SVGElement | null>>;
  STRANDS.forEach((strand) => {
    labelsStrandEnd[strand] = {} as Record<typeof ENDS[number], SVGElement | null>;
    ENDS.forEach((end) => {
      const isLabelActive = (strand === SENSE_STRAND) || (strand === ANTISENSE_STRAND && isAntisenseStrandActive);
      const labelText = STRAND_END_LABEL_TEXT[end][strand];
      const xPosition = xPositionOfStrandLabels[end];
      const labelY = POSITIONS_FOR_STRAND_ELEMENTS[strand].Y_POSITION_BASE_LABEL;

      labelsStrandEnd[strand][end] = isLabelActive
        ? svgElementFactory.createTextElement(labelText, xPosition, labelY, NUCLEOBASE_FONT_SIZE_PX, SVG_ELEMENT_COLORS.TEXT)
        : null;
    });
  });

  const labelsTerminusModification = {} as Record<typeof STRANDS[number], Record<typeof TERMINI[number], SVGElement | null>>;

  STRANDS.forEach((strand) => {
    labelsTerminusModification[strand] = {} as Record<typeof TERMINI[number], SVGElement | null>;

    TERMINI.forEach((terminus) => {
      if (strand === ANTISENSE_STRAND && !isAntisenseStrandActive) {
        labelsTerminusModification[strand][terminus] = null;
        return;
      }

      const end = (strand === SENSE_STRAND && terminus === FIVE_PRIME) ||
        (strand === ANTISENSE_STRAND && terminus === THREE_PRIME) ? END_LEFT : END_RIGHT;
      const labelY = POSITIONS_FOR_STRAND_ELEMENTS[strand].Y_POSITION_BASE_LABEL;
      const text = terminalModifications[strand][terminus];
      const xPosition = xPositionOfTerminusModifications[end][strand];

      labelsTerminusModification[strand][terminus] = svgElementFactory.createTextElement(
        text, xPosition, labelY, NUCLEOBASE_FONT_SIZE_PX, SVG_ELEMENT_COLORS.MODIFICATION_TEXT
      );
    });
  });

  const commentLabel = svgElementFactory.createTextElement(comment, xPositionOfStrandLabels[END_LEFT], SVG_VERTICAL_COORDINATES.COMMENT(isAntisenseStrandActive), COMMENT_FONT_SIZE_PX, SVG_ELEMENT_COLORS.TEXT);

  const starElement = isPhosphorothioateLinkageActive ? svgElementFactory.createStarElement(NUCLEOBASE_CIRCLE_RADIUS_PX, SVG_VERTICAL_COORDINATES.LEGEND_CIRCLES(isAntisenseStrandActive), SVG_ELEMENT_COLORS.LINKAGE_STAR) : null;

  const psLinkageLabel = isPhosphorothioateLinkageActive ? svgElementFactory.createTextElement('ps linkage', 2 * NUCLEOBASE_CIRCLE_RADIUS_PX - 8, SVG_VERTICAL_COORDINATES.LEGEND_TEXT(isAntisenseStrandActive), COMMENT_FONT_SIZE_PX, SVG_ELEMENT_COLORS.TEXT) : null;

  const svgElementsToAdd = [
    labelsStrandEnd[SENSE_STRAND][END_LEFT],
    labelsStrandEnd[SENSE_STRAND][END_RIGHT],
    labelsStrandEnd[ANTISENSE_STRAND][END_LEFT],
    labelsStrandEnd[ANTISENSE_STRAND][END_RIGHT],
    labelsTerminusModification[SENSE_STRAND][FIVE_PRIME],
    labelsTerminusModification[ANTISENSE_STRAND][THREE_PRIME],
    labelsTerminusModification[SENSE_STRAND][THREE_PRIME],
    labelsTerminusModification[ANTISENSE_STRAND][FIVE_PRIME],
    commentLabel,
    starElement,
    psLinkageLabel,
  ].filter((element) => element !== null) as SVGElement[];

  image.append(...svgElementsToAdd);

  const numberOfSsNucleotides = bases[SENSE_STRAND].filter((value) => !isOverhangNucleotide(value)).length;
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = bases[SENSE_STRAND].length - 1; i > -1; i--) {
    const xOfNumbers = computeNucleobaseCircleXPosition(i, rightOverhangs[SENSE_STRAND]) +
      computeNumericLabelOffset(bases[SENSE_STRAND], bases[SENSE_STRAND].length - i, numberOfSsNucleotides - nucleotideCounter);
    if (!isOverhangNucleotide(bases[SENSE_STRAND][i]))
      nucleotideCounter--;
    const nucleotideLabel = (!isOverhangNucleotide(bases[SENSE_STRAND][i]) && modificationsWithNumericLabels.includes(bases[SENSE_STRAND][i])) ?
      String(numberOfSsNucleotides - nucleotideCounter) : '';

    const baseCircleXCoord = computeNucleobaseCircleXPosition(i, rightOverhangs[SENSE_STRAND]);

    const nucleotideTextLabel = svgElementFactory.createTextElement(
      nucleotideLabel, xOfNumbers, POSITIONS_FOR_STRAND_ELEMENTS[SENSE_STRAND].Y_POSITION_NUMERIC_LABEL, COMMENT_FONT_SIZE_PX, SVG_ELEMENT_COLORS.TEXT
    );

    const nucleotideCircle = svgElementFactory.createCircleElement(
      baseCircleXCoord, POSITIONS_FOR_STRAND_ELEMENTS[SENSE_STRAND].Y_POSITION_NUCLEOBASE_CIRCLE, NUCLEOBASE_CIRCLE_RADIUS_PX, retrieveNucleobaseColor(bases[SENSE_STRAND][i])
    );

    const nucleotideBaseLabel = svgElementFactory.createTextElement(
      getNucleobaseLabelForCircle(bases[SENSE_STRAND], i), xOfNumbers, POSITIONS_FOR_STRAND_ELEMENTS[SENSE_STRAND].Y_POSITION_BASE_LABEL, NUCLEOBASE_FONT_SIZE_PX, determineTextColorBasedOnBackground(bases[SENSE_STRAND][i])
    );

    const linkageStar = phosphorothioateLinkages[SENSE_STRAND][i] ? svgElementFactory.createStarElement(
      baseCircleXCoord + NUCLEOBASE_CIRCLE_RADIUS_PX, POSITIONS_FOR_STRAND_ELEMENTS[SENSE_STRAND].Y_POSITION_BASE_LABEL + LINKAGE_STAR_RADIUS_PX, SVG_ELEMENT_COLORS.LINKAGE_STAR
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
      svgElementFactory.createStarElement(computeNucleobaseCircleXPosition(bases[SENSE_STRAND].length, rightOverhangs[SENSE_STRAND]) +
      NUCLEOBASE_CIRCLE_RADIUS_PX, POSITIONS_FOR_STRAND_ELEMENTS[SENSE_STRAND].Y_POSITION_BASE_LABEL + LINKAGE_STAR_RADIUS_PX, SVG_ELEMENT_COLORS.LINKAGE_STAR) : '',
  );

  const numberOfAsNucleotides = bases[ANTISENSE_STRAND].filter((value) => !isOverhangNucleotide(value)).length;
  if (isAntisenseStrandActive) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = bases[ANTISENSE_STRAND].length - 1; i > -1; i--) {
      if (!isOverhangNucleotide(bases[ANTISENSE_STRAND][i]))
        nucleotideCounter--;
      const xOfNumbers = computeNucleobaseCircleXPosition(i, rightOverhangs[ANTISENSE_STRAND]) +
        computeNumericLabelOffset(bases[ANTISENSE_STRAND], i, nucleotideCounter + 1);
      const n = (!isOverhangNucleotide(bases[ANTISENSE_STRAND][i]) && modificationsWithNumericLabels.includes(bases[ANTISENSE_STRAND][i])) ?
        String(nucleotideCounter + 1) : '';
      image.append(
        svgElementFactory.createTextElement(n, xOfNumbers, POSITIONS_FOR_STRAND_ELEMENTS[ANTISENSE_STRAND].Y_POSITION_NUMERIC_LABEL, COMMENT_FONT_SIZE_PX, SVG_ELEMENT_COLORS.TEXT),
        svgElementFactory.createCircleElement(computeNucleobaseCircleXPosition(i, rightOverhangs[ANTISENSE_STRAND]), POSITIONS_FOR_STRAND_ELEMENTS[ANTISENSE_STRAND].Y_POSITION_NUCLEOBASE_CIRCLE, NUCLEOBASE_CIRCLE_RADIUS_PX, retrieveNucleobaseColor(bases[ANTISENSE_STRAND][i])),
        svgElementFactory.createTextElement(getNucleobaseLabelForCircle(bases[ANTISENSE_STRAND], i),
          computeNucleobaseCircleXPosition(i, rightOverhangs[ANTISENSE_STRAND]) + computeNumericLabelOffset(bases[ANTISENSE_STRAND], i, nucleotideCounter + 1),
          POSITIONS_FOR_STRAND_ELEMENTS[ANTISENSE_STRAND].Y_POSITION_BASE_LABEL, NUCLEOBASE_FONT_SIZE_PX, determineTextColorBasedOnBackground(bases[ANTISENSE_STRAND][i])),
        phosphorothioateLinkages[ANTISENSE_STRAND][i] ? svgElementFactory.createStarElement(computeNucleobaseCircleXPosition(i, rightOverhangs[ANTISENSE_STRAND]) +
          NUCLEOBASE_CIRCLE_RADIUS_PX, POSITIONS_FOR_STRAND_ELEMENTS[ANTISENSE_STRAND].Y_POSITION_BASE_LABEL + LINKAGE_STAR_RADIUS_PX, SVG_ELEMENT_COLORS.LINKAGE_STAR) : '',
      );
    }
    image.append(
      phosphorothioateLinkages[ANTISENSE_STRAND][bases[ANTISENSE_STRAND].length] ?
        svgElementFactory.createStarElement(computeNucleobaseCircleXPosition(bases[ANTISENSE_STRAND].length, rightOverhangs[ANTISENSE_STRAND]) + NUCLEOBASE_CIRCLE_RADIUS_PX, POSITIONS_FOR_STRAND_ELEMENTS[ANTISENSE_STRAND].Y_POSITION_BASE_LABEL + LINKAGE_STAR_RADIUS_PX,
          SVG_ELEMENT_COLORS.LINKAGE_STAR) : '',
    );
  }

  const title = `${patternName} for ${numberOfSsNucleotides}${(isAntisenseStrandActive ? `/${numberOfAsNucleotides}` : '')}mer`;
  image.append(svgElementFactory.createTextElement(title, TITLE_POSITION.X, TITLE_POSITION.Y, NUCLEOBASE_FONT_SIZE_PX, SVG_ELEMENT_COLORS.TITLE_TEXT));
  for (let i = 0; i < distinctNucleobaseTypes.length; i++) {
    image.append(
      svgElementFactory.createCircleElement(computeLegendCircleXPosition(i), SVG_VERTICAL_COORDINATES.LEGEND_CIRCLES(isAntisenseStrandActive), LEGEND_CIRCLE_RADIUS_PX, retrieveNucleobaseColor(distinctNucleobaseTypes[i])),
      svgElementFactory.createTextElement(distinctNucleobaseTypes[i], computeLegendCircleXPosition(i) + LEGEND_CIRCLE_RADIUS_PX + 4, SVG_VERTICAL_COORDINATES.LEGEND_TEXT(isAntisenseStrandActive), COMMENT_FONT_SIZE_PX,
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
    return LIGHT_COLOR;
  }

  const [r, g, b] = rgbValues;
  const luminance = r * RED_COEFFICIENT + g * GREEN_COEFFICIENT + b * BLUE_COEFFICIENT;
  return luminance > LUMINANCE_THRESHOLD ? DARK_COLOR : LIGHT_COLOR;
}

function retrieveNucleobaseColor(nucleobase: string): string {
  return styleMap[nucleobase].color;
}
