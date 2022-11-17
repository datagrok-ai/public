import {NUCLEOTIDES} from '../structures-works/map';
import {isOverhang, svg, textWidth, countOverhangsOnTheRightEdge, getBaseColor, getTextInsideCircle,
  getFontColorVisibleOnBackground} from './helpers';

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
const modificationsColor = 'red';
const ssLeftText = 'SS: 5\'';
const asLeftText = 'AS: 3\'';
const ssRightText = '3\'';
const asRightText = '5\'';
const xOfTitle = BASE_RADIUS; // Math.round(width / 4),
const yOfTitle = BASE_RADIUS;
const xOfLeftTexts = 0;
const yOfSsNumbers = 2 * BASE_RADIUS;
const yOfAsNumbers = 8.5 * BASE_RADIUS;
const yOfSsTexts = 4 * BASE_RADIUS;
const yOfAsTexts = 7 * BASE_RADIUS;
const yOfSsCircles = 3.5 * BASE_RADIUS;
const yOfAsCircles = 6.5 * BASE_RADIUS;

export function drawAxolabsPattern(
  patternName: string, asExists: boolean, ssBases: string[],
  asBases: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[],
  ssThreeModification: string, ssFiveModification: string,
  asThreeModification: string, asFiveModification: string, comment: string,
  enumerateModifications: string[],
): Element {
  function getEquidistantXForLegend(index: number): number {
    return Math.round((index + startFrom) * width / (uniqueBases.length + startFrom) + LEGEND_RADIUS);
  }

  function getXOfBaseCircles(index: number, rightOverhangs: number): number {
    return widthOfRightModification +
      (resultingNumberOfNucleotidesInStrands - index + rightOverhangs + 1) * BASE_DIAMETER;
  }

  function getShiftToAlignNumberInsideCircle(bases: string[], generalIndex: number, nucleotideIndex: number): number {
    return (nucleotideIndex < 10 || NUCLEOTIDES.includes(bases[generalIndex])) ?
      shiftToAlignOneDigitNumberNearCircle : shiftToAlignTwoDigitNumberNearCircle;
  }

  ssBases = ssBases.reverse();
  ssPtoStatuses = ssPtoStatuses.reverse();

  const ssRightOverhangs = countOverhangsOnTheRightEdge(ssBases);
  const asRightOverhangs = countOverhangsOnTheRightEdge(asBases);
  const resultingNumberOfNucleotidesInStrands =
    Math.max(ssBases.length - ssRightOverhangs, asBases.length - asRightOverhangs);
  const widthOfBases =
    BASE_DIAMETER * (resultingNumberOfNucleotidesInStrands + Math.max(ssRightOverhangs, asRightOverhangs));
  const widthOfLeftModification =
    Math.max(textWidth(ssThreeModification, BASE_FONT_SIZE), textWidth(asFiveModification, BASE_FONT_SIZE));
  const widthOfRightModification =
    Math.max(textWidth(ssFiveModification, BASE_FONT_SIZE), textWidth(asThreeModification, BASE_FONT_SIZE));
  const widthOfLeftText = Math.max(textWidth(ssLeftText, BASE_FONT_SIZE), textWidth(asLeftText, BASE_FONT_SIZE));
  const widthOfRightText = Math.max(textWidth(ssRightText, BASE_FONT_SIZE), textWidth(asRightText, BASE_FONT_SIZE));
  const width =
    widthOfLeftText + widthOfLeftModification + widthOfBases +
    widthOfRightModification + widthOfRightText + BASE_DIAMETER;
  const height = asExists ? 11 * BASE_RADIUS : 9 * BASE_RADIUS;
  const uniqueBases = asExists ? [...new Set(ssBases.concat(asBases))] : [...new Set(ssBases)];
  const isPtoExist = asExists ?
    [...new Set(ssPtoStatuses.concat(asPtoStatuses))].includes(true) :
    [...new Set(ssPtoStatuses)].includes(true);
  const startFrom = isPtoExist ? 1 : 0;
  const xOfLeftModifications = xOfLeftTexts + widthOfLeftText - 5;
  const xOfSsRightModifications = ssRightOverhangs * BASE_DIAMETER + getXOfBaseCircles(-0.5, 0);
  const xOfAsRightModifications = asRightOverhangs * BASE_DIAMETER + getXOfBaseCircles(-0.5, 0);
  const xOfRightTexts =
    Math.max(xOfSsRightModifications, xOfAsRightModifications) +
      widthOfLeftModification + BASE_DIAMETER * (Math.max(ssRightOverhangs, asRightOverhangs));
  const yOfComment = asExists ? 11 * BASE_RADIUS : 8.5 * BASE_RADIUS;
  const yOfCirclesInLegends = asExists ? 9.5 * BASE_RADIUS : 6 * BASE_RADIUS;
  const yOfTextLegend = asExists ? 10 * BASE_RADIUS - 3 : yOfAsCircles - 3;

  const image = svg.render(width, height);

  image.append(
    svg.text(ssLeftText, xOfLeftTexts, yOfSsTexts, BASE_FONT_SIZE, FONT_COLOR),
    asExists ? svg.text(asLeftText, xOfLeftTexts, yOfAsTexts, BASE_FONT_SIZE, FONT_COLOR) : '',
    svg.text(ssRightText, xOfRightTexts, yOfSsTexts, BASE_FONT_SIZE, FONT_COLOR),
    asExists ? svg.text(asRightText, xOfRightTexts, yOfAsTexts, BASE_FONT_SIZE, FONT_COLOR) : '',
    svg.text(ssFiveModification, xOfLeftModifications, yOfSsTexts, BASE_FONT_SIZE, modificationsColor),
    asExists ? svg.text(asThreeModification, xOfLeftModifications, yOfAsTexts, BASE_FONT_SIZE, modificationsColor) : '',
    svg.text(ssThreeModification, xOfSsRightModifications, yOfSsTexts, BASE_FONT_SIZE, modificationsColor),
    asExists ? svg.text(asFiveModification, xOfAsRightModifications, yOfAsTexts, BASE_FONT_SIZE, modificationsColor) : '',
    svg.text(comment, xOfLeftTexts, yOfComment, LEGEND_FONT_SIZE, FONT_COLOR),
    isPtoExist ? svg.star(BASE_RADIUS, yOfCirclesInLegends, PS_LINKAGE_COLOR) : '',
    isPtoExist ? svg.text('ps linkage', 2 * BASE_RADIUS - 8, yOfTextLegend, LEGEND_FONT_SIZE, FONT_COLOR) : '',
  );

  let numberOfSsNucleotides = 0;
  for (let i = 0; i < ssBases.length; i++) {
    if (!isOverhang(ssBases[i]))
      numberOfSsNucleotides++;
  }
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = ssBases.length - 1; i > -1; i--) {
    const xOfNumbers = getXOfBaseCircles(i, ssRightOverhangs) +
      getShiftToAlignNumberInsideCircle(ssBases, ssBases.length - i, numberOfSsNucleotides - nucleotideCounter);
    if (!isOverhang(ssBases[i]))
      nucleotideCounter--;
    const n = (!isOverhang(ssBases[i]) && enumerateModifications.includes(ssBases[i])) ?
      String(numberOfSsNucleotides - nucleotideCounter) : '';
    image.append(
      svg.text(n, xOfNumbers, yOfSsNumbers, LEGEND_FONT_SIZE, FONT_COLOR),
      svg.circle(getXOfBaseCircles(i, ssRightOverhangs), yOfSsCircles, BASE_RADIUS, getBaseColor(ssBases[i])),
      svg.text(getTextInsideCircle(ssBases, i), xOfNumbers, yOfSsTexts, BASE_FONT_SIZE,
        getFontColorVisibleOnBackground(getBaseColor(ssBases[i]))),
      ssPtoStatuses[i] ?
        svg.star(getXOfBaseCircles(i, ssRightOverhangs) + BASE_RADIUS, yOfSsTexts + PS_LINKAGE_RADIUS, PS_LINKAGE_COLOR) :
        '',
    );
  }
  image.append(
    ssPtoStatuses[ssBases.length] ?
      svg.star(getXOfBaseCircles(ssBases.length, ssRightOverhangs) +
      BASE_RADIUS, yOfSsTexts + PS_LINKAGE_RADIUS, PS_LINKAGE_COLOR) : '',
  );

  let numberOfAsNucleotides = 0;
  for (let i = 0; i < asBases.length; i++) {
    if (!isOverhang(asBases[i]))
      numberOfAsNucleotides++;
  }
  if (asExists) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = asBases.length - 1; i > -1; i--) {
      if (!isOverhang(asBases[i]))
        nucleotideCounter--;
      const xOfNumbers = getXOfBaseCircles(i, asRightOverhangs) +
        getShiftToAlignNumberInsideCircle(asBases, i, nucleotideCounter + 1);
      const n = (!isOverhang(asBases[i]) && enumerateModifications.includes(asBases[i])) ?
        String(nucleotideCounter + 1) : '';
      image.append(
        svg.text(n, xOfNumbers, yOfAsNumbers, LEGEND_FONT_SIZE, FONT_COLOR),
        svg.circle(getXOfBaseCircles(i, asRightOverhangs), yOfAsCircles, BASE_RADIUS, getBaseColor(asBases[i])),
        svg.text(getTextInsideCircle(asBases, i),
          getXOfBaseCircles(i, asRightOverhangs) + getShiftToAlignNumberInsideCircle(asBases, i, nucleotideCounter + 1),
          yOfAsTexts, BASE_FONT_SIZE, getFontColorVisibleOnBackground(getBaseColor(asBases[i]))),
        asPtoStatuses[i] ? svg.star(getXOfBaseCircles(i, asRightOverhangs) +
          BASE_RADIUS, yOfAsTexts + PS_LINKAGE_RADIUS, PS_LINKAGE_COLOR) : '',
      );
    }
    image.append(
      asPtoStatuses[asBases.length] ?
        svg.star(getXOfBaseCircles(asBases.length, asRightOverhangs) +
         BASE_RADIUS, yOfAsTexts + PS_LINKAGE_RADIUS, PS_LINKAGE_COLOR) : '',
    );
  }

  const title = `${patternName} for ${numberOfSsNucleotides}${(asExists ? `/${numberOfAsNucleotides}` : '')}mer`;
  image.append(svg.text(title, xOfTitle, yOfTitle, BASE_FONT_SIZE, TITLE_FONT_COLOR));
  for (let i = 0; i < uniqueBases.length; i++) {
    image.append(
      svg.circle(getEquidistantXForLegend(i), yOfCirclesInLegends, LEGEND_RADIUS, getBaseColor(uniqueBases[i])),
      svg.text(uniqueBases[i], getEquidistantXForLegend(i) +
        LEGEND_RADIUS + 4, yOfTextLegend, LEGEND_FONT_SIZE, FONT_COLOR),
    );
  }
  return image;
}
