import { study } from "../clinical-study";
import * as ui from "datagrok-api/ui";
import { domains } from "../sdtm-meta";

export function updateDivInnerHTML(div: HTMLDivElement, content: any){
    div.innerHTML = '';
    div.append(content);
  }

export function checkMissingDomains(requiredDomains: string[], or: boolean, obj: any) {
  let missingDomains = requiredDomains.filter(it => study.domains[it] === null);
  let noMissingDomains = or ? requiredDomains.some(it => study.domains[it] !== null) : !missingDomains.length; 
  if(noMissingDomains){
    obj.createView();
  } else {
    let domainsDiv = ui.div();
    missingDomains.forEach(it => {domainsDiv.append(ui.divText(it))})
    obj.root.append(ui.div([
      ui.h2('Missing domains:'),
      domainsDiv
    ], {style: {margin: 'auto', textAlign: 'center'}}));
  }
}