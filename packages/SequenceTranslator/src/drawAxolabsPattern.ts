/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {axolabsMap} from "./axolabsMap";

export function drawAxolabsPattern(patternName: string, ssBaseStatuses: string[], asBaseStatuses: string[], ssPtoStatuses: boolean[], asPtoStatuses: boolean[]) {
  const svg = {
    xmlns : "http://www.w3.org/2000/svg",
    render : function(width: string, height: string) {
      const e = document.createElementNS(this.xmlns, 'svg');
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
    radius = 20,
    width = (2 * radius + 1) * maxNumberInStrands,
    height = 200,
    title = patternName + ' for ' + String(ssBaseStatuses.length) + '/' + String(asBaseStatuses.length) + 'mer',
    fontSize = '30',
    fontWeight = 'normal',
    fontColor = 'black';

  let image = svg.render(String(width), String(height));
  image.append(svg.text(title, String(Math.round(width / 2)), String(2 * radius), fontSize, fontWeight, fontColor));
  image.append(svg.text('SS:', '0', String(4 * radius), fontSize, fontWeight, fontColor));
  image.append(svg.text('AS:', '0', String(7 * radius), fontSize, fontWeight, fontColor));
  for (let i = 0; i < ssBaseStatuses.length; i++) {
    image.append(
      svg.circle(String((i + 2) * 2 * radius), String(4 * radius), String(radius), axolabsMap[ssBaseStatuses[i]]["color"])
    );
  }
  for (let i = 0; i < asBaseStatuses.length; i++) {
    image.append(
      svg.circle(String((i + 2) * 2 * radius), String(7 * radius), String(radius), axolabsMap[asBaseStatuses[i]]["color"])
    );
  }
  return image;
}