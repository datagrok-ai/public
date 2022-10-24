import {axolabsMap} from './constants';

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
    points += currentX + ',' + currentY + ' ';
  }
  return points;
}

function countOverhangsOnTheRightEdge(modifications: string[]): number {
  let i = 0;
  while (i < modifications.length && modifications[i].slice(-3) == '(o)')
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

function getTextInsideCircle(bases: string[], index: number, enumerateModifications: string[]): string {
  return (
    bases[index].slice(-3) == '(o)' ||
    !enumerateModifications.includes(bases[index]) ||
    !['A', 'G', 'C', 'U', 'T'].includes(bases[index])
  ) ?
    '' :
    bases[index];
}

function getFontColorVisibleOnBackground(rgbString: string): string {
  const rgbIntList = rgbString.match(/\d+/g)!.map((e) => Number(e));
  return (rgbIntList[0] * 0.299 + rgbIntList[1] * 0.587 + rgbIntList[2] * 0.114) > 186 ? '#33333' : '#ffffff';
}

function getBaseColor(base: string): string {
  return axolabsMap[base]['color'];
}

export function drawAxolabsPattern(
  patternName: string, asExists: boolean, ssBases: string[],
  asBases: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[],
  ssThreeModification: string, ssFiveModification: string,
  asThreeModification: string, asFiveModification: string, comment: string,
  enumerateModifications: string[]) {
  function getEquidistantXForLegend(index: number): number {
    return Math.round((index + startFrom) * width / (uniqueBases.length + startFrom) + legendRadius);
  }

  function getXOfBaseCircles(index: number, rightOverhangs: number): number {
    return widthOfRightModification +
      (resultingNumberOfNucleotidesInStrands - index + rightOverhangs + 1) * baseDiameter;
  }

  function getShiftToAlignNumberInsideCircle(bases: string[], generalIndex: number, nucleotideIndex: number): number {
    return (nucleotideIndex < 10 || ['A', 'G', 'C', 'U', 'T'].
      includes(bases[generalIndex])) ? shiftToAlignOneDigitNumberInsideCircle : shiftToAlignTwoDigitNumberInsideCircle;
  }

  const svg = {
    xmlns: 'http://www.w3.org/2000/svg',
    render: function(width: number, height: number) {
      const e = document.createElementNS(this.xmlns, 'svg');
      e.setAttribute('id', 'mySvg');
      e.setAttribute('width', String(width));
      e.setAttribute('height', String(height));
      return e;
    },
    circle: function(x: number, y: number, radius: number, color: string) {
      const e = document.createElementNS(this.xmlns, 'circle');
      e.setAttribute('cx', String(x));
      e.setAttribute('cy', String(y));
      e.setAttribute('r', String(radius));
      e.setAttribute('fill', color);
      return e;
    },
    text: function(text: string, x: number, y: number, fontSize: number, color: string) {
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
    star: function(x: number, y: number, fill: string) {
      const e = document.createElementNS(this.xmlns, 'polygon');
      e.setAttribute('points', getPointsToDrawStar(x, y));
      e.setAttribute('fill', fill);
      return e;
    },
  };

  ssBases = ssBases.reverse();
  ssPtoStatuses = ssPtoStatuses.reverse();

  const baseRadius = 15;
  const shiftToAlignTwoDigitNumberInsideCircle = -10;
  const shiftToAlignOneDigitNumberInsideCircle = -5;
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

  const ssRightOverhangs = countOverhangsOnTheRightEdge(ssBases);
  const asRightOverhangs = countOverhangsOnTheRightEdge(asBases);
  const resultingNumberOfNucleotidesInStrands =
    Math.max(ssBases.length - ssRightOverhangs, asBases.length - asRightOverhangs);
  const baseDiameter = 2 * baseRadius;
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
  const height = asExists ? 11 * baseRadius : 9 * baseRadius;
  const xOfTitle = baseRadius; // Math.round(width / 4),
  const uniqueBases = asExists ? [...new Set(ssBases.concat(asBases))] : [...new Set(ssBases)];
  const isPtoExist =
    asExists ? [...new Set(ssPtoStatuses.concat(asPtoStatuses))].includes(true) :
      [...new Set(ssPtoStatuses)].includes(true);
  const startFrom = isPtoExist ? 1 : 0;
  const xOfLeftTexts = 0;
  const xOfLeftModifications = xOfLeftTexts + widthOfLeftText - 5;
  const xOfSsRightModifications = ssRightOverhangs * baseDiameter + getXOfBaseCircles(-0.5, 0);
  const xOfAsRightModifications = asRightOverhangs * baseDiameter + getXOfBaseCircles(-0.5, 0);
  const xOfRightTexts =
    Math.max(xOfSsRightModifications, xOfAsRightModifications) +
      widthOfLeftModification + baseDiameter * (Math.max(ssRightOverhangs, asRightOverhangs));
  const yOfTitle = baseRadius;
  const yOfSsNumbers = 2 * baseRadius;
  const yOfAsNumbers = 8.5 * baseRadius;
  const yOfSsTexts = 4 * baseRadius;
  const yOfAsTexts = 7 * baseRadius;
  const yOfComment = asExists ? 11 * baseRadius : 8.5 * baseRadius;
  const yOfSsCircles = 3.5 * baseRadius;
  const yOfAsCircles = 6.5 * baseRadius;
  const yOfCirclesInLegends = asExists ? 9.5 * baseRadius : 6 * baseRadius;
  const yOfTextLegend = asExists ? 10 * baseRadius - 3 : yOfAsCircles - 3;

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
    isPtoExist ? svg.star(baseRadius, yOfCirclesInLegends, psLinkageColor) : '',
    isPtoExist ? svg.text('ps linkage', 2 * baseRadius - 8, yOfTextLegend, legendFontSize, fontColor) : '',
  );

  let numberOfSsNucleotides = 0;
  for (let i = 0; i < ssBases.length; i++) {
    if (ssBases[i].slice(-3) != '(o)')
      numberOfSsNucleotides++;
  }
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = ssBases.length - 1; i > -1; i--) {
    const xOfNumbers = getXOfBaseCircles(i, ssRightOverhangs) +
      getShiftToAlignNumberInsideCircle(ssBases, ssBases.length - i, numberOfSsNucleotides - nucleotideCounter);
    if (ssBases[i].slice(-3) != '(o)')
      nucleotideCounter--;
    image.append(
      svg.text(String(numberOfSsNucleotides - nucleotideCounter), xOfNumbers, yOfSsNumbers, legendFontSize, fontColor),
      svg.circle(getXOfBaseCircles(i, ssRightOverhangs), yOfSsCircles, baseRadius, getBaseColor(ssBases[i])),
      svg.text(getTextInsideCircle(ssBases, i, enumerateModifications),
        xOfNumbers, yOfSsTexts, baseFontSize, getFontColorVisibleOnBackground(getBaseColor(ssBases[i]))),
      ssPtoStatuses[i] ?
        svg.star(getXOfBaseCircles(i, ssRightOverhangs) + baseRadius, yOfSsTexts + psLinkageRadius, psLinkageColor) :
        '',
    );
  }
  image.append(
    ssPtoStatuses[ssBases.length] ?
      svg.star(getXOfBaseCircles(ssBases.length, ssRightOverhangs) +
      baseRadius, yOfSsTexts + psLinkageRadius, psLinkageColor) : '',
  );

  let numberOfAsNucleotides = 0;
  for (let i = 0; i < asBases.length; i++) {
    if (asBases[i].slice(-3) != '(o)')
      numberOfAsNucleotides++;
  }
  if (asExists) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = asBases.length - 1; i > -1; i--) {
      if (asBases[i].slice(-3) != '(o)')
        nucleotideCounter--;
      const xOfNumbers = getXOfBaseCircles(i, asRightOverhangs) +
        getShiftToAlignNumberInsideCircle(asBases, asBases.length - i, numberOfAsNucleotides - nucleotideCounter);
      image.append(
        svg.text(String(numberOfAsNucleotides - nucleotideCounter), xOfNumbers, yOfAsNumbers, legendFontSize, fontColor),
        svg.circle(getXOfBaseCircles(i, asRightOverhangs), yOfAsCircles, baseRadius, getBaseColor(asBases[i])),
        svg.text(getTextInsideCircle(asBases, i, enumerateModifications),
          getXOfBaseCircles(i, asRightOverhangs) + getShiftToAlignNumberInsideCircle(asBases, i, nucleotideCounter + 1),
          yOfAsTexts, baseFontSize, getFontColorVisibleOnBackground(getBaseColor(asBases[i]))),
        asPtoStatuses[i] ? svg.star(getXOfBaseCircles(i, asRightOverhangs) +
          baseRadius, yOfAsTexts + psLinkageRadius, psLinkageColor) : '',
      );
    }
    image.append(
      asPtoStatuses[asBases.length] ?
        svg.star(getXOfBaseCircles(asBases.length, asRightOverhangs) +
         baseRadius, yOfAsTexts + psLinkageRadius, psLinkageColor) : '',
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
