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

function getTextWidth(text: string, font: number): number {
  const context = document.createElement('canvas').getContext('2d');
  // @ts-ignore
  context.font = String(font);
  // @ts-ignore
  return 2 * context.measureText(text).width;
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

  function getXOfBaseCircles(index: number): number {
    return widthOfRightModification + (maxNumberOfNucleotidesInStrands - index + 1) * baseDiameter;
  }

  function getShiftToAlignNumberInsideCircle(index: number): number {
    return index < 10 ? shiftToAlignOneDigitNumberInsideCircle : shiftToAlignTwoDigitNumberInsideCircle;
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

  const maxNumberOfNucleotidesInStrands = Math.max(ssBases.length, asBases.length),
    baseDiameter = 2 * baseRadius,
    widthOfBases = baseDiameter * maxNumberOfNucleotidesInStrands,
    widthOfLeftModification = Math.max(getTextWidth(ssThreeModification, baseFontSize), getTextWidth(asFiveModification, baseFontSize)),
    widthOfRightModification = Math.max(getTextWidth(ssFiveModification, baseFontSize), getTextWidth(asThreeModification, baseFontSize)),
    widthOfLeftText = Math.max(getTextWidth(ssLeftText, baseFontSize), getTextWidth(asLeftText, baseFontSize)),
    widthOfRightText = Math.max(getTextWidth(ssRightText, baseFontSize), getTextWidth(asRightText, baseFontSize)),
    width = widthOfLeftText + widthOfLeftModification + widthOfBases + widthOfRightModification + widthOfRightText + baseDiameter,
    height = asExists ? 11 * baseRadius : 9 * baseRadius,
    title = patternName + ' for ' + String(ssBases.length) + (asExists ? '/' + String(asBases.length) : '') + 'mer',
    xOfTitle = Math.round(width / 2),
    uniqueBases = asExists ? [...new Set(ssBases.concat(asBases))] : [...new Set(ssBases)],
    isPtoExist = asExists ? [...new Set(ssPtoStatuses.concat(asPtoStatuses))].includes(true) : [...new Set(ssPtoStatuses)].includes(true),
    startFrom = isPtoExist ? 1 : 0,
    xOfLeftTexts = 0,
    xOfLeftModifications = xOfLeftTexts + widthOfLeftText - 5,
    xOfRightModifications = getXOfBaseCircles(-0.5),
    xOfRightTexts = xOfRightModifications + widthOfLeftModification,
    yOfTitle = 2 * baseRadius,
    yOfSsTexts = 4 * baseRadius,
    yOfAsTexts = 7 * baseRadius,
    xOfComment = 2 * baseRadius,
    yOfComment = asExists ? 10.5 * baseRadius : 8 * baseRadius,
    yOfSsCircles = 3.5 * baseRadius,
    yOfAsCircles = 6.5 * baseRadius,
    yOfCirclesInLegends = asExists ? 9 * baseRadius : 6 * baseRadius,
    yOfTextLegend = asExists ? 9.5 * baseRadius - 3 : yOfAsCircles - 3;

  let image = svg.render(width, height);

  image.append(
    svg.text(title, xOfTitle, yOfTitle, baseFontSize, fontColor),
    svg.text(ssLeftText, xOfLeftTexts, yOfSsTexts, baseFontSize, fontColor),
    asExists ? svg.text(asLeftText, xOfLeftTexts, yOfAsTexts, baseFontSize, fontColor) : '',
    svg.text(ssRightText, xOfRightTexts, yOfSsTexts, baseFontSize, fontColor),
    asExists ? svg.text(asRightText, xOfRightTexts, yOfAsTexts, baseFontSize, fontColor) : '',
    svg.text(ssFiveModification, xOfLeftModifications, yOfSsTexts, baseFontSize, modificationsColor),
    asExists ? svg.text(asThreeModification, xOfLeftModifications, yOfAsTexts, baseFontSize, modificationsColor) : '',
    svg.text(ssThreeModification, xOfRightModifications, yOfSsTexts, baseFontSize, modificationsColor),
    asExists ? svg.text(asFiveModification, xOfRightModifications, yOfAsTexts, baseFontSize, modificationsColor) : '',
    svg.text(comment, xOfComment, yOfComment, baseFontSize, fontColor),
    isPtoExist ? svg.star(baseRadius, yOfCirclesInLegends, psLinkageColor) : '',
    isPtoExist ? svg.text('ps linkage', 2 * baseRadius - 8, yOfTextLegend, legendFontSize, fontColor) : ''
  );

  for (let i = ssBases.length - 1; i > -1; i--)
    image.append(
      svg.circle(getXOfBaseCircles(i), yOfSsCircles, baseRadius, getBaseColor(ssBases[i])),
      svg.text(String(ssBases.length - i), getXOfBaseCircles(i) + getShiftToAlignNumberInsideCircle(ssBases.length - i), yOfSsTexts, baseFontSize, getFontColorVisibleOnBackground(axolabsMap[ssBases[i]]["color"])),
      ssPtoStatuses[i] ? svg.star(getXOfBaseCircles(i) + baseRadius, yOfSsTexts + psLinkageRadius, psLinkageColor) : ''
    );
  image.append(
    ssPtoStatuses[ssBases.length] ? svg.star(getXOfBaseCircles(ssBases.length) + baseRadius, yOfSsTexts + psLinkageRadius, psLinkageColor) : ''
  );

  if (asExists) {
    for (let i = asBases.length - 1; i > -1; i--)
      image.append(
        svg.circle(getXOfBaseCircles(i), yOfAsCircles, baseRadius, getBaseColor(asBases[i])),
        svg.text(String(i + 1), getXOfBaseCircles(i) + getShiftToAlignNumberInsideCircle(i + 1), yOfAsTexts, baseFontSize, getFontColorVisibleOnBackground(axolabsMap[asBases[i]]["color"])),
        asPtoStatuses[i] ? svg.star(getXOfBaseCircles(i) + baseRadius, yOfAsTexts + psLinkageRadius, psLinkageColor) : ''
      );
    image.append(
      asPtoStatuses[asBases.length] ? svg.star(getXOfBaseCircles(asBases.length) + baseRadius, yOfAsTexts + psLinkageRadius, psLinkageColor) : ''
    );
  }

  for (let i = 0; i < uniqueBases.length; i++)
    image.append(
      svg.circle(getEquidistantXForLegend(i), yOfCirclesInLegends, legendRadius, getBaseColor(uniqueBases[i])),
      svg.text(uniqueBases[i], getEquidistantXForLegend(i) + legendRadius, yOfTextLegend, legendFontSize, fontColor)
    );
  return image;
}