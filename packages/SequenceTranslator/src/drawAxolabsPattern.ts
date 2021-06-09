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

export function drawAxolabsPattern(patternName: string, createAsStrand: boolean, ssBaseStatuses: string[], asBaseStatuses: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[], ssModif: string, asModif:string) {
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
    polygon : function (x: number, y: number, fill: string) {
      const e = document.createElementNS(this.xmlns,"polygon");
      e.setAttribute("points", getStarPoints(x, y));
      e.setAttribute('fill', fill);
      return e;
    }
  };

  const maxNumberInStrands = Math.max(ssBaseStatuses.length, asBaseStatuses.length),
    baseRadius = 15,
    psLinkageRadius = 5,
    psLinkageColor = 'red',
    fontSize = '17',
    width = (2 * baseRadius + 1) * maxNumberInStrands + 5 * baseRadius + Math.max(ssModif.length, asModif.length) * 9,
    height = 11 * baseRadius,
    title = patternName + ' for ' + String(ssBaseStatuses.length) + '/' + String(asBaseStatuses.length) + 'mer',
    textShift = 5,
    fontWeight = 'normal',
    fontColor = 'grey';

  let image = svg.render(String(width), String(height));
  image.append(svg.text(title, String(Math.round(width / 2)), String(2 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text('SS:', '0', String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text('AS:', '0', String(7 * baseRadius), fontSize, fontWeight, fontColor));

  image.append(svg.text("3'", String(2.5 * baseRadius), String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("5'", String(width - baseRadius), String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("3'", String(width - baseRadius), String(7 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("5'", String(2.5 * baseRadius), String(7 * baseRadius), fontSize, fontWeight, fontColor));

  image.append(svg.text(ssModif, String((maxNumberInStrands + 2) * 2 * baseRadius + baseRadius), String(4 * baseRadius), fontSize, fontWeight, psLinkageColor));
  image.append(svg.text(asModif, String((maxNumberInStrands + 2) * 2 * baseRadius + baseRadius), String(7 * baseRadius), fontSize, fontWeight, psLinkageColor));

  for (let i = ssBaseStatuses.length - 1; i > -1; i--) {
    image.append(
      svg.circle(String((maxNumberInStrands - i + 2) * 2 * baseRadius), String(3.5 * baseRadius), String(baseRadius), axolabsMap[ssBaseStatuses[i]]["color"])
    );
    image.append(
      svg.text(String(ssBaseStatuses.length - i), String((maxNumberInStrands - i + 2) * 2 * baseRadius - baseRadius + textShift), String(4 * baseRadius), fontSize, fontWeight, fontColor)
    );
    if (ssPtoStatuses[i]) {
      image.append(
        svg.polygon((maxNumberInStrands - i + 2) * 2 * baseRadius + baseRadius, 4 * baseRadius + psLinkageRadius, psLinkageColor)
      );
    }
  }
  if (createAsStrand) {
    for (let i = asBaseStatuses.length - 1; i > -1; i--) {
      image.append(
        svg.circle(String((maxNumberInStrands - i + 2) * 2 * baseRadius), String(6.5 * baseRadius), String(baseRadius), axolabsMap[asBaseStatuses[i]]["color"])
      );
      image.append(
        svg.text(String(i + 1), String((maxNumberInStrands - i + 2) * 2 * baseRadius - baseRadius + textShift), String(7 * baseRadius), fontSize, fontWeight, fontColor)
      )
      if (asPtoStatuses[i]) {
        image.append(
          svg.polygon((maxNumberInStrands - i + 2) * 2 * baseRadius - baseRadius, 7 * baseRadius + psLinkageRadius, psLinkageColor)
        );
      }
    }
  }
  const uniqueBases = [...new Set(ssBaseStatuses.concat(asBaseStatuses))];
  for (let i = 0; i < uniqueBases.length; i++) {
    image.append(
      svg.circle(String(Math.round(i * width / uniqueBases.length + baseRadius)), String(9.5 * baseRadius), String(baseRadius), axolabsMap[uniqueBases[i]]["color"])
    );
    image.append(
      svg.text(uniqueBases[i], String(Math.round(i * width / uniqueBases.length) + 2 * baseRadius), String(10 * baseRadius), fontSize, fontWeight, fontColor)
    );
  }
  return image;
}