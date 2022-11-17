import {AXOLABS_MAP} from './constants';
import {NUCLEOTIDES} from '../structures-works/map';
import {isOverhang} from './helpers';

// https://uxdesign.cc/star-rating-make-svg-great-again-d4ce4731347e
function getPointsToDrawStar(centerX: number, centerY: number): string {
  const innerCirclePoints = 5; // a 5 point star
  const innerRadius = 15 / innerCirclePoints;
  const innerOuterRadiusRatio = 2; // outter circle is x2 the inner
  const outerRadius = innerRadius * innerOuterRadiusRatio;
  const angle = Math.PI / innerCirclePoints;
  const angleOffsetToCenterStar = 60;
  const totalNumberOfPoints = innerCirclePoints * 2; // 10 in a 5-points star

  let points = '';
  for (let i = 0; i < totalNumberOfPoints; i++) {
    const r = (i % 2 == 0) ? outerRadius : innerRadius;
    const currentX = centerX + Math.cos(i * angle + angleOffsetToCenterStar) * r;
    const currentY = centerY + Math.sin(i * angle + angleOffsetToCenterStar) * r;
    points += `${currentX},${currentY} `;
  }
  return points;
}

function countOverhangsOnTheRightEdge(modifications: string[]): number {
  let i = 0;
  while (i < modifications.length && isOverhang(modifications[i]))
    i++;
  return (i == modifications.length - 1) ? 0 : i;
}

function getTextWidth(text: string, font: number): number {
  const context = document.createElement('canvas').getContext('2d');
  // @ts-ignore
  context.font = String(font);
  // @ts-ignore
  return 2 * context.measureText(text).width;
}

function getTextInsideCircle(bases: string[], index: number): string {
  return (isOverhang(bases[index]) || !NUCLEOTIDES.includes(bases[index])) ? '' : bases[index];
}

function getFontColorVisibleOnBackground(rgbString: string): string {
  const rgbIntList = rgbString.match(/\d+/g)!.map((e) => Number(e));
  return (rgbIntList[0] * 0.299 + rgbIntList[1] * 0.587 + rgbIntList[2] * 0.114) > 186 ? '#33333' : '#ffffff';
}

function getBaseColor(base: string): string {
  return AXOLABS_MAP[base]['color'];
}

const svg = {
  xmlns: 'http://www.w3.org/2000/svg',
  render: function(width: number, height: number): Element {
    const e = document.createElementNS(this.xmlns, 'svg');
    e.setAttribute('id', 'mySvg');
    e.setAttribute('width', String(width));
    e.setAttribute('height', String(height));
    return e;
  },
  circle: function(x: number, y: number, radius: number, color: string): Element {
    const e = document.createElementNS(this.xmlns, 'circle');
    e.setAttribute('cx', String(x));
    e.setAttribute('cy', String(y));
    e.setAttribute('r', String(radius));
    e.setAttribute('fill', color);
    return e;
  },
  text: function(text: string, x: number, y: number, fontSize: number, color: string): Element {
    const e = document.createElementNS(this.xmlns, 'text');
    e.setAttribute('x', String(x));
    e.setAttribute('y', String(y));
    e.setAttribute('font-size', String(fontSize));
    e.setAttribute('font-weight', 'normal');
    e.setAttribute('font-family', 'Arial');
    e.setAttribute('fill', color);
    e.innerHTML = text;
    return e;
  },
  star: function(x: number, y: number, fill: string): Element {
    const e = document.createElementNS(this.xmlns, 'polygon');
    e.setAttribute('points', getPointsToDrawStar(x, y));
    e.setAttribute('fill', fill);
    return e;
  },
};

const BASE_RADIUS = 15;
const shiftToAlignTwoDigitNumberNearCircle = -10;
const shiftToAlignOneDigitNumberNearCircle = -5;
const legendRadius = 6;
const psLinkageRadius = 5;
const baseFontSize = 17;
const legendFontSize = 14;
const psLinkageColor = 'red';
const fontColor = 'var(--grey-6)';
const titleFontColor = 'black';
const modificationsColor = 'red';
const ssLeftText = 'SS: 5\'';
const asLeftText = 'AS: 3\'';
const ssRightText = '3\'';
const asRightText = '5\'';

export function drawAxolabsPattern(
  patternName: string, asExists: boolean, ssBases: string[],
  asBases: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[],
  ssThreeModification: string, ssFiveModification: string,
  asThreeModification: string, asFiveModification: string, comment: string,
  enumerateModifications: string[],
): Element {
  function getEquidistantXForLegend(index: number): number {
    return Math.round((index + startFrom) * width / (uniqueBases.length + startFrom) + legendRadius);
  }

  function getXOfBaseCircles(index: number, rightOverhangs: number): number {
    return widthOfRightModification +
      (resultingNumberOfNucleotidesInStrands - index + rightOverhangs + 1) * baseDiameter;
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
  const baseDiameter = 2 * BASE_RADIUS;
  const widthOfBases =
    baseDiameter * (resultingNumberOfNucleotidesInStrands + Math.max(ssRightOverhangs, asRightOverhangs));
  const widthOfLeftModification =
    Math.max(getTextWidth(ssThreeModification, baseFontSize), getTextWidth(asFiveModification, baseFontSize));
  const widthOfRightModification =
    Math.max(getTextWidth(ssFiveModification, baseFontSize), getTextWidth(asThreeModification, baseFontSize));
  const widthOfLeftText = Math.max(getTextWidth(ssLeftText, baseFontSize), getTextWidth(asLeftText, baseFontSize));
  const widthOfRightText = Math.max(getTextWidth(ssRightText, baseFontSize), getTextWidth(asRightText, baseFontSize));
  const width =
    widthOfLeftText + widthOfLeftModification + widthOfBases +
    widthOfRightModification + widthOfRightText + baseDiameter;
  const height = asExists ? 11 * BASE_RADIUS : 9 * BASE_RADIUS;
  const xOfTitle = BASE_RADIUS; // Math.round(width / 4),
  const uniqueBases = asExists ? [...new Set(ssBases.concat(asBases))] : [...new Set(ssBases)];
  const isPtoExist = asExists ?
    [...new Set(ssPtoStatuses.concat(asPtoStatuses))].includes(true) :
    [...new Set(ssPtoStatuses)].includes(true);
  const startFrom = isPtoExist ? 1 : 0;
  const xOfLeftTexts = 0;
  const xOfLeftModifications = xOfLeftTexts + widthOfLeftText - 5;
  const xOfSsRightModifications = ssRightOverhangs * baseDiameter + getXOfBaseCircles(-0.5, 0);
  const xOfAsRightModifications = asRightOverhangs * baseDiameter + getXOfBaseCircles(-0.5, 0);
  const xOfRightTexts =
    Math.max(xOfSsRightModifications, xOfAsRightModifications) +
      widthOfLeftModification + baseDiameter * (Math.max(ssRightOverhangs, asRightOverhangs));
  const yOfTitle = BASE_RADIUS;
  const yOfSsNumbers = 2 * BASE_RADIUS;
  const yOfAsNumbers = 8.5 * BASE_RADIUS;
  const yOfSsTexts = 4 * BASE_RADIUS;
  const yOfAsTexts = 7 * BASE_RADIUS;
  const yOfComment = asExists ? 11 * BASE_RADIUS : 8.5 * BASE_RADIUS;
  const yOfSsCircles = 3.5 * BASE_RADIUS;
  const yOfAsCircles = 6.5 * BASE_RADIUS;
  const yOfCirclesInLegends = asExists ? 9.5 * BASE_RADIUS : 6 * BASE_RADIUS;
  const yOfTextLegend = asExists ? 10 * BASE_RADIUS - 3 : yOfAsCircles - 3;

  const image = svg.render(width, height);

  image.append(
    svg.text(ssLeftText, xOfLeftTexts, yOfSsTexts, baseFontSize, fontColor),
    asExists ? svg.text(asLeftText, xOfLeftTexts, yOfAsTexts, baseFontSize, fontColor) : '',
    svg.text(ssRightText, xOfRightTexts, yOfSsTexts, baseFontSize, fontColor),
    asExists ? svg.text(asRightText, xOfRightTexts, yOfAsTexts, baseFontSize, fontColor) : '',
    svg.text(ssFiveModification, xOfLeftModifications, yOfSsTexts, baseFontSize, modificationsColor),
    asExists ? svg.text(asThreeModification, xOfLeftModifications, yOfAsTexts, baseFontSize, modificationsColor) : '',
    svg.text(ssThreeModification, xOfSsRightModifications, yOfSsTexts, baseFontSize, modificationsColor),
    asExists ? svg.text(asFiveModification, xOfAsRightModifications, yOfAsTexts, baseFontSize, modificationsColor) : '',
    svg.text(comment, xOfLeftTexts, yOfComment, legendFontSize, fontColor),
    isPtoExist ? svg.star(BASE_RADIUS, yOfCirclesInLegends, psLinkageColor) : '',
    isPtoExist ? svg.text('ps linkage', 2 * BASE_RADIUS - 8, yOfTextLegend, legendFontSize, fontColor) : '',
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
      svg.text(n, xOfNumbers, yOfSsNumbers, legendFontSize, fontColor),
      svg.circle(getXOfBaseCircles(i, ssRightOverhangs), yOfSsCircles, BASE_RADIUS, getBaseColor(ssBases[i])),
      svg.text(getTextInsideCircle(ssBases, i), xOfNumbers, yOfSsTexts, baseFontSize,
        getFontColorVisibleOnBackground(getBaseColor(ssBases[i]))),
      ssPtoStatuses[i] ?
        svg.star(getXOfBaseCircles(i, ssRightOverhangs) + BASE_RADIUS, yOfSsTexts + psLinkageRadius, psLinkageColor) :
        '',
    );
  }
  image.append(
    ssPtoStatuses[ssBases.length] ?
      svg.star(getXOfBaseCircles(ssBases.length, ssRightOverhangs) +
      BASE_RADIUS, yOfSsTexts + psLinkageRadius, psLinkageColor) : '',
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
        svg.text(n, xOfNumbers, yOfAsNumbers, legendFontSize, fontColor),
        svg.circle(getXOfBaseCircles(i, asRightOverhangs), yOfAsCircles, BASE_RADIUS, getBaseColor(asBases[i])),
        svg.text(getTextInsideCircle(asBases, i),
          getXOfBaseCircles(i, asRightOverhangs) + getShiftToAlignNumberInsideCircle(asBases, i, nucleotideCounter + 1),
          yOfAsTexts, baseFontSize, getFontColorVisibleOnBackground(getBaseColor(asBases[i]))),
        asPtoStatuses[i] ? svg.star(getXOfBaseCircles(i, asRightOverhangs) +
          BASE_RADIUS, yOfAsTexts + psLinkageRadius, psLinkageColor) : '',
      );
    }
    image.append(
      asPtoStatuses[asBases.length] ?
        svg.star(getXOfBaseCircles(asBases.length, asRightOverhangs) +
         BASE_RADIUS, yOfAsTexts + psLinkageRadius, psLinkageColor) : '',
    );
  }

  const title = patternName + ' for ' +
    String(numberOfSsNucleotides) + (asExists ? '/' + String(numberOfAsNucleotides) : '') + 'mer';
  image.append(svg.text(title, xOfTitle, yOfTitle, baseFontSize, titleFontColor));
  for (let i = 0; i < uniqueBases.length; i++) {
    image.append(
      svg.circle(getEquidistantXForLegend(i), yOfCirclesInLegends, legendRadius, getBaseColor(uniqueBases[i])),
      svg.text(uniqueBases[i], getEquidistantXForLegend(i) +
        legendRadius + 4, yOfTextLegend, legendFontSize, fontColor),
    );
  }
  return image;
}
