/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {axolabsMap} from "./axolabsMap";

export function drawAxolabsPattern(patternName: string, createAsStrand: boolean, ssBaseStatuses: string[], asBaseStatuses: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[]) {
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
    }
  };

  const maxNumberInStrands = Math.max(ssBaseStatuses.length, asBaseStatuses.length),
    baseRadius = 20,
    psLinkageRadius = 5,
    psLinkageColor = 'red',
    width = (2 * baseRadius + 1) * maxNumberInStrands + 5 * baseRadius,
    height = 11 * baseRadius,
    title = patternName + ' for ' + String(ssBaseStatuses.length) + '/' + String(asBaseStatuses.length) + 'mer',
    fontSize = '30',
    textShift = 3,
    fontWeight = 'normal',
    fontColor = 'black';

  let image = svg.render(String(width), String(height));
  image.append(svg.text(title, String(Math.round(width / 2)), String(2 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text('SS:', '0', String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text('AS:', '0', String(7 * baseRadius), fontSize, fontWeight, fontColor));

  image.append(svg.text("3'", String(2.5 * baseRadius), String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("5'", String(width - baseRadius), String(4 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("3'", String(width - baseRadius), String(7 * baseRadius), fontSize, fontWeight, fontColor));
  image.append(svg.text("5'", String(2.5 * baseRadius), String(7 * baseRadius), fontSize, fontWeight, fontColor));

  for (let i = ssBaseStatuses.length - 1; i > -1; i--) {
    image.append(
      svg.circle(String((maxNumberInStrands - i + 2) * 2 * baseRadius), String(3.5 * baseRadius), String(baseRadius), axolabsMap[ssBaseStatuses[i]]["color"])
    );
    image.append(
      svg.text(String(ssBaseStatuses.length - i), String((maxNumberInStrands - i + 2) * 2 * baseRadius - baseRadius + textShift), String(4 * baseRadius), fontSize, fontWeight, fontColor)
    );
    if (ssPtoStatuses[i]) {
      image.append(
        svg.circle(String((maxNumberInStrands - i + 2) * 2 * baseRadius + baseRadius), String(4 * baseRadius + psLinkageRadius), String(psLinkageRadius), psLinkageColor)
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
          svg.circle(String((maxNumberInStrands - i + 2) * 2 * baseRadius + baseRadius), String(7 * baseRadius + psLinkageRadius), String(psLinkageRadius), psLinkageColor)
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