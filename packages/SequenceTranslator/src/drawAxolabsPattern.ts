import {axolabsMap} from "./axolabsMap";

// https://uxdesign.cc/star-rating-make-svg-great-again-d4ce4731347e
function getPointsToDrawStar(centerX: number, centerY: number) {
  const innerCirclePoints = 5;  // a 5 point star
  const innerRadius = 15 / innerCirclePoints;
  const innerOuterRadiusRatio = 2; // outter circle is x2 the inner
  const outerRadius = innerRadius * innerOuterRadiusRatio;
  const angle = Math.PI / innerCirclePoints;
  const angleOffsetToCenterStar = 60;
  const totalNumberOfPoints = innerCirclePoints * 2; // 10 in a 5-points star

  let points = '';
  for (let i = 0; i < totalNumberOfPoints; i++) {
    let r = (i % 2 == 0) ? outerRadius : innerRadius;
    let currentX = centerX + Math.cos(i * angle + angleOffsetToCenterStar) * r;
    let currentY = centerY + Math.sin(i * angle + angleOffsetToCenterStar) * r;
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

function getTextInsideCircle(bases: string[], index: number, nucleotideCounter: number, numberOfNucleotides: number): string {
  return (bases[index].slice(-3) == "(o)") ? "" :
    ['A', 'G', 'C', 'U', 'T'].includes(bases[index]) ? bases[index] : String(numberOfNucleotides - nucleotideCounter);
}

function getFontColorVisibleOnBackground(rgbString: string) {
  const rgbIntList = rgbString.match(/\d+/g)!.map(e => Number(e));
  return (rgbIntList[0] * 0.299 + rgbIntList[1] * 0.587 + rgbIntList[2] * 0.114) > 186 ? '#333333' : '#ffffff';
}

function getBaseColor(base: string): string {
  return axolabsMap[base]["color"];
}

export function drawAxolabsPattern(patternName: string, asExists: boolean, ssBases: string[], asBases: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[],
                                   ssThreeModification: string, ssFiveModification: string, asThreeModification: string, asFiveModification: string, comment: string) {

  function getEquidistantXForLegend(index: number): number {
    return Math.round((index + startFrom) * width / (uniqueBases.length + startFrom) + legendRadius);
  }

  function getXOfBaseCircles(index: number, rightOverhangs: number): number {
    return widthOfRightModification + (resultingNumberOfNucleotidesInStrands - index + rightOverhangs + 1) * baseDiameter;
  }

  function getShiftToAlignNumberInsideCircle(bases: string[], generalIndex: number, nucleotideIndex: number): number {
    return (nucleotideIndex < 10 || ['A', 'G', 'C', 'U', 'T'].includes(bases[generalIndex])) ? shiftToAlignOneDigitNumberInsideCircle : shiftToAlignTwoDigitNumberInsideCircle;
  }

  const svg = {
    xmlns: "http://www.w3.org/2000/svg",
    render: function(width: number, height: number) {
      const e = document.createElementNS(this.xmlns, 'svg');
      e.setAttribute("id", "mySvg");
      e.setAttribute("width", String(width));
      e.setAttribute("height", String(height));
      return e;
    },
    circle: function(x: number, y: number, radius: number, color: string) {
      const e = document.createElementNS(this.xmlns, 'circle');
      e.setAttribute("cx", String(x));
      e.setAttribute("cy", String(y));
      e.setAttribute("r", String(radius));
      e.setAttribute("fill", color);
      return e;
    },
    text: function(text: string, x: number, y: number, fontSize: number, color: string) {
      const e = document.createElementNS(this.xmlns, 'text');
      e.setAttribute("x", String(x));
      e.setAttribute("y", String(y));
      e.setAttribute("font-size", String(fontSize));
      e.setAttribute("font-weight", "normal");
      e.setAttribute("fill", color);
      e.innerHTML = text;
      return e;
    },
    star: function(x: number, y: number, fill: string) {
      const e = document.createElementNS(this.xmlns, "polygon");
      e.setAttribute("points", getPointsToDrawStar(x, y));
      e.setAttribute("fill", fill);
      return e;
    }
  };

  ssBases = ssBases.reverse();
  ssPtoStatuses = ssPtoStatuses.reverse();

  const baseRadius = 15,
    shiftToAlignTwoDigitNumberInsideCircle = -10,
    shiftToAlignOneDigitNumberInsideCircle = -5,
    legendRadius = 6,
    psLinkageRadius = 5,
    baseFontSize = 17,
    legendFontSize = 14,
    psLinkageColor = 'red',
    fontColor = 'var(--grey-6)',
    modificationsColor = 'red',
    ssLeftText = "SS: 5'",
    asLeftText = "AS: 3'",
    ssRightText = "3'",
    asRightText = "5'";

  const ssRightOverhangs = countOverhangsOnTheRightEdge(ssBases);
  const asRightOverhangs = countOverhangsOnTheRightEdge(asBases);
  const resultingNumberOfNucleotidesInStrands = Math.max(ssBases.length - ssRightOverhangs, asBases.length - asRightOverhangs),
    baseDiameter = 2 * baseRadius,
    widthOfBases = baseDiameter * (resultingNumberOfNucleotidesInStrands + Math.max(ssRightOverhangs, asRightOverhangs)),
    widthOfLeftModification = Math.max(getTextWidth(ssThreeModification, baseFontSize), getTextWidth(asFiveModification, baseFontSize)),
    widthOfRightModification = Math.max(getTextWidth(ssFiveModification, baseFontSize), getTextWidth(asThreeModification, baseFontSize)),
    widthOfLeftText = Math.max(getTextWidth(ssLeftText, baseFontSize), getTextWidth(asLeftText, baseFontSize)),
    widthOfRightText = Math.max(getTextWidth(ssRightText, baseFontSize), getTextWidth(asRightText, baseFontSize)),
    width = widthOfLeftText + widthOfLeftModification + widthOfBases + widthOfRightModification + widthOfRightText + baseDiameter,
    height = asExists ? 11 * baseRadius : 9 * baseRadius,
    xOfTitle = Math.round(width / 4),
    uniqueBases = asExists ? [...new Set(ssBases.concat(asBases))] : [...new Set(ssBases)],
    isPtoExist = asExists ? [...new Set(ssPtoStatuses.concat(asPtoStatuses))].includes(true) : [...new Set(ssPtoStatuses)].includes(true),
    startFrom = isPtoExist ? 1 : 0,
    xOfLeftTexts = 0,
    xOfLeftModifications = xOfLeftTexts + widthOfLeftText - 5,
    xOfSsRightModifications = ssRightOverhangs * baseDiameter + getXOfBaseCircles(-0.5, 0),
    xOfAsRightModifications = asRightOverhangs * baseDiameter + getXOfBaseCircles(-0.5, 0),
    xOfRightTexts = Math.max(xOfSsRightModifications, xOfAsRightModifications) + widthOfLeftModification + baseDiameter * (Math.max(ssRightOverhangs, asRightOverhangs)),
    yOfTitle = 2 * baseRadius,
    yOfSsTexts = 4 * baseRadius,
    yOfAsTexts = 7 * baseRadius,
    yOfComment = asExists ? 11 * baseRadius : 8.5 * baseRadius,
    yOfSsCircles = 3.5 * baseRadius,
    yOfAsCircles = 6.5 * baseRadius,
    yOfCirclesInLegends = asExists ? 9 * baseRadius : 6 * baseRadius,
    yOfTextLegend = asExists ? 9.5 * baseRadius - 3 : yOfAsCircles - 3;

  let image = svg.render(width, height);

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
    isPtoExist ? svg.text('ps linkage', 2 * baseRadius - 8, yOfTextLegend, legendFontSize, fontColor) : ''
  );

  let numberOfSsNucleotides = 0;
  for (let i = 0; i < ssBases.length; i++)
    if (ssBases[i].slice(-3) != '(o)')
      numberOfSsNucleotides++;
  let nucleotideCounter = numberOfSsNucleotides;
  for (let i = ssBases.length - 1; i > -1; i--) {
    if (ssBases[i].slice(-3) != '(o)')
      nucleotideCounter--;
    image.append(
      svg.circle(getXOfBaseCircles(i, ssRightOverhangs), yOfSsCircles, baseRadius, getBaseColor(ssBases[i])),
      svg.text(getTextInsideCircle(ssBases, i, nucleotideCounter, numberOfSsNucleotides), getXOfBaseCircles(i, ssRightOverhangs) + getShiftToAlignNumberInsideCircle(ssBases, ssBases.length - i, numberOfSsNucleotides - nucleotideCounter), yOfSsTexts, baseFontSize, getFontColorVisibleOnBackground(axolabsMap[ssBases[i]]["color"])),
      ssPtoStatuses[i] ? svg.star(getXOfBaseCircles(i, ssRightOverhangs) + baseRadius, yOfSsTexts + psLinkageRadius, psLinkageColor) : ''
    );
  }
  image.append(
    ssPtoStatuses[ssBases.length] ? svg.star(getXOfBaseCircles(ssBases.length, ssRightOverhangs) + baseRadius, yOfSsTexts + psLinkageRadius, psLinkageColor) : ''
  );

  let numberOfAsNucleotides = 0;
  for (let i = 0; i < asBases.length; i++)
    if (asBases[i].slice(-3) != '(o)')
      numberOfAsNucleotides++;
  if (asExists) {
    let nucleotideCounter = numberOfAsNucleotides;
    for (let i = asBases.length - 1; i > -1; i--) {
      if (asBases[i].slice(-3) != '(o)')
        nucleotideCounter--;
      image.append(
        svg.circle(getXOfBaseCircles(i, asRightOverhangs), yOfAsCircles, baseRadius, getBaseColor(asBases[i])),
        svg.text(getTextInsideCircle(asBases, i, numberOfAsNucleotides - nucleotideCounter - 1, numberOfAsNucleotides), getXOfBaseCircles(i, asRightOverhangs) + getShiftToAlignNumberInsideCircle(asBases, i, nucleotideCounter + 1), yOfAsTexts, baseFontSize, getFontColorVisibleOnBackground(axolabsMap[asBases[i]]["color"])),
        asPtoStatuses[i] ? svg.star(getXOfBaseCircles(i, asRightOverhangs) + baseRadius, yOfAsTexts + psLinkageRadius, psLinkageColor) : ''
      );
    }
    image.append(
      asPtoStatuses[asBases.length] ? svg.star(getXOfBaseCircles(asBases.length, asRightOverhangs) + baseRadius, yOfAsTexts + psLinkageRadius, psLinkageColor) : ''
    );
  }

  let title = patternName + ' for ' + String(numberOfSsNucleotides) + (asExists ? '/' + String(numberOfAsNucleotides) : '') + 'mer';
  image.append(svg.text(title, xOfTitle, yOfTitle, baseFontSize, fontColor));
  for (let i = 0; i < uniqueBases.length; i++)
    image.append(
      svg.circle(getEquidistantXForLegend(i), yOfCirclesInLegends, legendRadius, getBaseColor(uniqueBases[i])),
      svg.text(uniqueBases[i], getEquidistantXForLegend(i) + legendRadius, yOfTextLegend, legendFontSize, fontColor)
    );
  return image;
}