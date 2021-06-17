/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {axolabsMap} from "./axolabsMap";

function getStarPoints(centerX: number, centerY: number) {
  const innerCirclePoints = 5;  // a 5 point star
  const innerRadius = 15 / innerCirclePoints;
  const innerOuterRadiusRatio = 2; // outter circle is x2 the inner
  const outerRadius = innerRadius * innerOuterRadiusRatio;
  const angle = (Math.PI / innerCirclePoints);
  const angleOffsetToCenterStar = 60;
  const totalPoints = innerCirclePoints * 2; // 10 in a 5-points star

  let points = '';
  for (let i = 0; i < totalPoints; i++) {
    let isEvenIndex = i % 2 == 0;
    let r = isEvenIndex ? outerRadius : innerRadius;
    let currX = centerX + Math.cos(i * angle + angleOffsetToCenterStar ) * r;
    let currY = centerY + Math.sin(i * angle + angleOffsetToCenterStar) * r;
    points += currX + ',' + currY + ' ';
  }
  return points;
}

function getFontColor(rgbString: string) {
  let rgbIntList = rgbString.match(/\d+/g)!.map(e => Number(e));
  return ((rgbIntList[0] * 0.299 + rgbIntList[1] * 0.587 + rgbIntList[2] * 0.114) > 186) ? '#333333' : '#ffffff';
}

export function drawAxolabsPattern(patternName: string, createAsStrand: boolean, ssBaseStatuses: string[], asBaseStatuses: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[],
                                   sSthreeModification: string, sSfiveModification: string, aSthreeModification: string, aSfiveModification: string, comment: string) {
  ssBaseStatuses = ssBaseStatuses.reverse();
  ssPtoStatuses = ssPtoStatuses.reverse();

  const svg = {
    xmlns : "http://www.w3.org/2000/svg",
    render : function(width: string, height: string) {
      const e = document.createElementNS(this.xmlns, 'svg');
      e.setAttribute('id', 'mySvg');
      e.setAttribute("width", width);
      e.setAttribute("height", height);
      return e;
    },
    circle : function(x: string, y: string, r: string, fill: string) {
      const e = document.createElementNS(this.xmlns, 'circle');
      e.setAttribute("cx", x);
      e.setAttribute("cy", y);
      e.setAttribute("r", r);
      e.setAttribute("fill", fill);
      return e;
    },
    text : function(text: string, x: string, y: string, fontSize: string, fontWeight: string, fill: string) {
      const e = document.createElementNS(this.xmlns, 'text');
      e.setAttribute("x", x);
      e.setAttribute("y", y);
      e.setAttribute("font-size", fontSize);
      e.setAttribute("font-weight", fontWeight);
      e.setAttribute("fill", fill);
      e.innerHTML = text;
      return e;
    },
    polygon : function(x: number, y: number, fill: string) {
      const e = document.createElementNS(this.xmlns, "polygon");
      e.setAttribute("points", getStarPoints(x, y));
      e.setAttribute('fill', fill);
      return e;
    }
  };

  const maxNumberInStrands = Math.max(ssBaseStatuses.length, asBaseStatuses.length),
    baseRadius = 15,
    legendRadius = 6,
    legendFontSize = '14',
    psLinkageRadius = 5,
    psLinkageColor = 'red',
    fontSize = '17',
    width = (2 * baseRadius + 1) * maxNumberInStrands + 5 * baseRadius + (Math.max(sSthreeModification.length, aSfiveModification.length) + Math.max(sSfiveModification.length, aSthreeModification.length)) * 10,
    height = 11 * baseRadius,
    title = patternName + ' for ' + String(ssBaseStatuses.length) + ((createAsStrand) ? '/' + String(asBaseStatuses.length) : '') + 'mer',
    textShift = 5,
    fontWeight = 'normal',
    fontColor = 'var(--grey-6)',
    shift = Math.max(sSfiveModification.length, aSthreeModification.length) * 10;

  let image = svg.render(String(width), String(height));
  image.append(svg.text(title, String(Math.round(width / 2)), String(2 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text('SS:', '0', String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text('AS:', '0', String(7 * baseRadius), fontSize, fontWeight, fontColor));

  image.append(svg.text("5'", String(2.5 * baseRadius), String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("3'", String(width - 2 * legendRadius), String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("5'", String(width - 2 * legendRadius), String(7 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("3'", String(2.5 * baseRadius), String(7 * baseRadius), fontSize, fontWeight, fontColor));

  image.append(svg.text(sSfiveModification, String(Math.round(3.5 * baseRadius)), String(4 * baseRadius), fontSize, fontWeight, psLinkageColor));
  image.append(svg.text(aSthreeModification, String(Math.round(3.5 * baseRadius)), String(7 * baseRadius), fontSize, fontWeight, psLinkageColor));

  image.append(svg.text(sSthreeModification, String(shift + (maxNumberInStrands + 2.5) * 2 * baseRadius), String(4 * baseRadius), fontSize, fontWeight, psLinkageColor));
  image.append(svg.text(aSfiveModification, String(shift + (maxNumberInStrands + 2.5) * 2 * baseRadius), String(7 * baseRadius), fontSize, fontWeight, psLinkageColor));

  for (let i = ssBaseStatuses.length - 1; i > -1; i--) {
    image.append(
      svg.circle(String(shift + (maxNumberInStrands - i + 2) * 2 * baseRadius), String(3.5 * baseRadius), String(baseRadius), axolabsMap[ssBaseStatuses[i]]["color"])
    );
    image.append(
      svg.text(String(ssBaseStatuses.length - i), String(shift + (maxNumberInStrands - i + 2) * 2 * baseRadius - baseRadius + ((ssBaseStatuses.length - i < 10) ? textShift + 5 : textShift)), String(4 * baseRadius), fontSize, fontWeight, getFontColor(axolabsMap[ssBaseStatuses[i]]["color"]))
    );
    if (ssPtoStatuses[i]) {
      image.append(
        svg.polygon(shift + (maxNumberInStrands - i + 2) * 2 * baseRadius + baseRadius, 4 * baseRadius + psLinkageRadius, psLinkageColor)
      );
    }
  }
  if (createAsStrand) {
    for (let i = asBaseStatuses.length - 1; i > -1; i--) {
      image.append(
        svg.circle(String(shift + (maxNumberInStrands - i + 2) * 2 * baseRadius), String(6.5 * baseRadius), String(baseRadius), axolabsMap[asBaseStatuses[i]]["color"])
      );
      image.append(
        svg.text(String(i + 1), String(shift + (maxNumberInStrands - i + 2) * 2 * baseRadius - baseRadius + ((i < 9) ? textShift + 5 : textShift)), String(7 * baseRadius), fontSize, fontWeight, getFontColor(axolabsMap[asBaseStatuses[i]]["color"]))
      )
      if (asPtoStatuses[i]) {
        image.append(
          svg.polygon(shift + (maxNumberInStrands - i + 2) * 2 * baseRadius - baseRadius, 7 * baseRadius + psLinkageRadius, psLinkageColor)
        );
      }
    }
  }
  const uniqueBases = [...new Set(ssBaseStatuses.concat(asBaseStatuses))];
  const isPtoExist = [...new Set(ssPtoStatuses.concat(asPtoStatuses))].includes(true);
  const startFrom = isPtoExist ? 1 : 0;
  if (isPtoExist) {
    image.append(
      svg.polygon(baseRadius, 9 * baseRadius, psLinkageColor)
    );
    image.append(
      svg.text('ps linkage', String(2 * baseRadius - 8), String(9.5 * baseRadius - 3), legendFontSize, fontWeight, fontColor)
    );
  }
  for (let i = 0; i < uniqueBases.length; i++) {
    image.append(
      svg.circle(String(Math.round((i + startFrom) * width / (uniqueBases.length + startFrom) + baseRadius)), String(9 * baseRadius), String(legendRadius), axolabsMap[uniqueBases[i]]["color"])
    );
    image.append(
      svg.text(uniqueBases[i], String(Math.round((i + startFrom) * width / (uniqueBases.length + startFrom)) + 2 * baseRadius - 8), String(9.5 * baseRadius - 3), legendFontSize, fontWeight, fontColor)
    );
  }
  image.append(svg.text(comment, String(2 * baseRadius), String(10.5 * baseRadius), fontSize, fontWeight, fontColor));
  return image;
}