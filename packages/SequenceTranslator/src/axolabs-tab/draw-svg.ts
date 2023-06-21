import {NUCLEOTIDES} from '../hardcode-to-be-eliminated/map';
import {isOverhang, svg, textWidth, countOverhangsOnTheRightEdge, baseColor, textInsideCircle,
  fontColorVisibleOnBackground, isOneDigitNumber} from './helpers';

const BASE_RADIUS = 15;
const BASE_DIAMETER = 2 * BASE_RADIUS;
const shiftToAlignTwoDigitNumberNearCircle = -10;
const shiftToAlignOneDigitNumberNearCircle = -5;
const LEGEND_RADIUS = 6;
const PS_LINKAGE_RADIUS = 5;
const BASE_FONT_SIZE = 17;
const LEGEND_FONT_SIZE = 14;
const PS_LINKAGE_COLOR = 'red';
const FONT_COLOR = 'var(--grey-6)';
const TITLE_FONT_COLOR = 'black';
const MODIFICATIONS_COLOR = 'red';
const SS_LEFT_TEXT = 'SS: 5\'';
const AS_LEFT_TEXT = 'AS: 3\'';
const SS_RIGHT_TEXT = '3\'';
const AS_RIGHT_TEXT = '5\'';

const WIDTH_OF_LEFT_TEXT = Math.max(
  textWidth(SS_LEFT_TEXT, BASE_FONT_SIZE),
  textWidth(AS_LEFT_TEXT, BASE_FONT_SIZE),
);

const WIDTH_OF_RIGHT_TEXT = Math.max(
  textWidth(SS_RIGHT_TEXT, BASE_FONT_SIZE),
  textWidth(AS_RIGHT_TEXT, BASE_FONT_SIZE),
);

const X = {
  TITLE: BASE_RADIUS, // Math.round(width / 4),
  LEFT_TEXTS: 0,
};
const X_OF_LEFT_MODIFICATIONS = X.LEFT_TEXTS + WIDTH_OF_LEFT_TEXT - 5;

const Y = {
  TITLE: BASE_RADIUS,
  SS_INDICES: 2 * BASE_RADIUS,
  SS_CIRCLES: 3.5 * BASE_RADIUS,
  SS_TEXTS: 4 * BASE_RADIUS,
  AS_CIRCLES: 6.5 * BASE_RADIUS,
  AS_TEXTS: 7 * BASE_RADIUS,
  AS_INDICES: 8.5 * BASE_RADIUS,
  comment: (asExists: boolean) => (asExists) ? 11 * BASE_RADIUS : 8.5 * BASE_RADIUS,
  circlesInLegends: (asExists: boolean) => (asExists) ? 9.5 * BASE_RADIUS : 6 * BASE_RADIUS,
  textLegend: (asExists: boolean) => (asExists) ? 10 * BASE_RADIUS - 3 : Y.AS_CIRCLES - 3,
  svgHeight: (asExists: boolean) => (asExists) ? 11 * BASE_RADIUS : 9 * BASE_RADIUS,
};

export function drawAxolabsPattern(
  patternName: string, asExists: boolean, ssBases: string[],
  asBases: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[],
  ss3Modification: string, ss5Modification: string,
  as3Modification: string, as5Modification: string, comment: string,
  enumerateModifications: string[],
): Element {
  function equidistantXForLegend(index: number): number {
    return Math.round((index + startFrom) * width / (uniqueBases.length + startFrom) + LEGEND_RADIUS);
  }

  function xOfBaseCircles(index: number, rightOverhangs: number): number {
    return widthOfRightModification +
      (resultingNumberOfNucleotidesInStrands - index + rightOverhangs + 1) * BASE_DIAMETER;
  }

  function shiftToAlignNumberNearCircle(bases: string[], generalIndex: number, nucleotideIndex: number): number {
    return (isOneDigitNumber(nucleotideIndex) || NUCLEOTIDES.includes(bases[generalIndex])) ?
      shiftToAlignOneDigitNumberNearCircle : shiftToAlignTwoDigitNumberNearCircle;
  }

  ssBases = ssBases.reverse();
  ssPtoStatuses = ssPtoStatuses.reverse();

  const ssRightOverhangs = countOverhangsOnTheRightEdge(ssBases);
  const asRightOverhangs = countOverhangsOnTheRightEdge(asBases);

  const resultingNumberOfNucleotidesInStrands = Math.max(
    ssBases.length - ssRightOverhangs,
    asBases.length - asRightOverhangs,
  );

  const widthOfRightOverhangs = Math.max(ssRightOverhangs, asRightOverhangs);
  const widthOfBases = BASE_DIAMETER * (resultingNumberOfNucleotidesInStrands + widthOfRightOverhangs);

  const widthOfLeftModification = Math.max(
    textWidth(ss3Modification, BASE_FONT_SIZE),
    textWidth(as5Modification, BASE_FONT_SIZE),
  );

  const widthOfRightModification = Math.max(
    textWidth(ss5Modification, BASE_FONT_SIZE),
    textWidth(as3Modification, BASE_FONT_SIZE),
  );

  const uniqueBases = asExists ?
    [...new Set(ssBases.concat(asBases))] :
    [...new Set(ssBases)];

  const isPtoExist = asExists ?
    ssPtoStatuses.concat(asPtoStatuses).includes(true) :
    ssPtoStatuses.includes(true);

  const startFrom = isPtoExist ? 1 : 0;

  const xOfSsRightModifications = ssRightOverhangs * BASE_DIAMETER + xOfBaseCircles(-0.5, 0);
  const xOfAsRightModifications = asRightOverhangs * BASE_DIAMETER + xOfBaseCircles(-0.5, 0);

  const xOfRightTexts = Math.max(xOfSsRightModifications, xOfAsRightModifications) + widthOfLeftModification +
    BASE_DIAMETER * widthOfRightOverhangs;

  const width = WIDTH_OF_LEFT_TEXT + widthOfLeftModification + widthOfBases + widthOfRightModification +
                WIDTH_OF_RIGHT_TEXT + BASE_DIAMETER;
  const image = svg.render(width, Y.svgHeight(asExists));

  image.append(
    svg.text(SS_LEFT_TEXT, X.LEFT_TEXTS, Y.SS_TEXTS, BASE_FONT_SIZE, FONT_COLOR),
    asExists ? svg.text(AS_LEFT_TEXT, X.LEFT_TEXTS, Y.AS_TEXTS, BASE_FONT_SIZE, FONT_COLOR) : '',
    svg.text(SS_RIGHT_TEXT, xOfRightTexts, Y.SS_TEXTS, BASE_FONT_SIZE, FONT_COLOR),
    asExists ? svg.text(AS_RIGHT_TEXT, xOfRightTexts, Y.AS_TEXTS, BASE_FONT_SIZE, FONT_COLOR) : '',
    svg.text(ss5Modification, X_OF_LEFT_MODIFICATIONS, Y.SS_TEXTS, BASE_FONT_SIZE, MODIFICATIONS_COLOR),
    asExists ? svg.text(as3Modification, X_OF_LEFT_MODIFICATIONS, Y.AS_TEXTS, BASE_FONT_SIZE, MODIFICATIONS_COLOR) : '',
    svg.text(ss3Modification, xOfSsRightModifications, Y.SS_TEXTS, BASE_FONT_SIZE, MODIFICATIONS_COLOR),
    asExists ? svg.text(as5Modification, xOfAsRightModifications, Y.AS_TEXTS, BASE_FONT_SIZE, MODIFICATIONS_COLOR) : '',
    svg.text(comment, X.LEFT_TEXTS, Y.comment(asExists), LEGEND_FONT_SIZE, FONT_COLOR),
    isPtoExist ? svg.star(BASE_RADIUS, Y.circlesInLegends(asExists), PS_LINKAGE_COLOR) : '',
    isPtoExist ? svg.text('ps linkage', 2 * BASE_RADIUS - 8, Y.textLegend(asExists), LEGEND_FONT_SIZE, FONT_COLOR) : '',
  );

  const numberOfSsNucleotides = ssBases.filter((value) => !isOverhang(value)).length;
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = ssBases.length - 1; i > -1; i--) {
    const xOfNumbers = xOfBaseCircles(i, ssRightOverhangs) +
      shiftToAlignNumberNearCircle(ssBases, ssBases.length - i, numberOfSsNucleotides - nucleotideCounter);
    if (!isOverhang(ssBases[i]))
      nucleotideCounter--;
    const n = (!isOverhang(ssBases[i]) && enumerateModifications.includes(ssBases[i])) ?
      String(numberOfSsNucleotides - nucleotideCounter) : '';
    image.append(
      svg.text(n, xOfNumbers, Y.SS_INDICES, LEGEND_FONT_SIZE, FONT_COLOR),
      svg.circle(xOfBaseCircles(i, ssRightOverhangs), Y.SS_CIRCLES, BASE_RADIUS, baseColor(ssBases[i])),
      svg.text(textInsideCircle(ssBases, i), xOfNumbers, Y.SS_TEXTS, BASE_FONT_SIZE,
        fontColorVisibleOnBackground(ssBases[i])),
      ssPtoStatuses[i] ?
        svg.star(xOfBaseCircles(i, ssRightOverhangs) + BASE_RADIUS, Y.SS_TEXTS + PS_LINKAGE_RADIUS, PS_LINKAGE_COLOR) :
        '',
    );
  }
  image.append(
    ssPtoStatuses[ssBases.length] ?
      svg.star(xOfBaseCircles(ssBases.length, ssRightOverhangs) +
      BASE_RADIUS, Y.SS_TEXTS + PS_LINKAGE_RADIUS, PS_LINKAGE_COLOR) : '',
  );

  const numberOfAsNucleotides = asBases.filter((value) => !isOverhang(value)).length;
  if (asExists) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = asBases.length - 1; i > -1; i--) {
      if (!isOverhang(asBases[i]))
        nucleotideCounter--;
      const xOfNumbers = xOfBaseCircles(i, asRightOverhangs) +
        shiftToAlignNumberNearCircle(asBases, i, nucleotideCounter + 1);
      const n = (!isOverhang(asBases[i]) && enumerateModifications.includes(asBases[i])) ?
        String(nucleotideCounter + 1) : '';
      image.append(
        svg.text(n, xOfNumbers, Y.AS_INDICES, LEGEND_FONT_SIZE, FONT_COLOR),
        svg.circle(xOfBaseCircles(i, asRightOverhangs), Y.AS_CIRCLES, BASE_RADIUS, baseColor(asBases[i])),
        svg.text(textInsideCircle(asBases, i),
          xOfBaseCircles(i, asRightOverhangs) + shiftToAlignNumberNearCircle(asBases, i, nucleotideCounter + 1),
          Y.AS_TEXTS, BASE_FONT_SIZE, fontColorVisibleOnBackground(asBases[i])),
        asPtoStatuses[i] ? svg.star(xOfBaseCircles(i, asRightOverhangs) +
          BASE_RADIUS, Y.AS_TEXTS + PS_LINKAGE_RADIUS, PS_LINKAGE_COLOR) : '',
      );
    }
    image.append(
      asPtoStatuses[asBases.length] ?
        svg.star(xOfBaseCircles(asBases.length, asRightOverhangs) + BASE_RADIUS, Y.AS_TEXTS + PS_LINKAGE_RADIUS,
          PS_LINKAGE_COLOR) : '',
    );
  }

  const title = `${patternName} for ${numberOfSsNucleotides}${(asExists ? `/${numberOfAsNucleotides}` : '')}mer`;
  image.append(svg.text(title, X.TITLE, Y.TITLE, BASE_FONT_SIZE, TITLE_FONT_COLOR));
  for (let i = 0; i < uniqueBases.length; i++) {
    image.append(
      svg.circle(equidistantXForLegend(i), Y.circlesInLegends(asExists), LEGEND_RADIUS, baseColor(uniqueBases[i])),
      svg.text(uniqueBases[i], equidistantXForLegend(i) + LEGEND_RADIUS + 4, Y.textLegend(asExists), LEGEND_FONT_SIZE,
        FONT_COLOR),
    );
  }
  return image;
}
