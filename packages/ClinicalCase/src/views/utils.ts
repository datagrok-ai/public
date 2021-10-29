import { study } from "../clinical-study";
import * as ui from "datagrok-api/ui";

export function updateDivInnerHTML(div: HTMLDivElement, content: any){
    div.innerHTML = '';
    div.append(content);
  }

export function checkDomainExists(requiredDomains: string[], or: boolean, obj: any) {
  let domainsExist = or ? requiredDomains.some(it => study.domains[it] !== null) : requiredDomains.every(it => study.domains[it] !== null)
  if(domainsExist){
    obj.createView();
  } else {
    let oneOfStr = or ? 'At least one of' : 'Each of';
    requiredDomains = requiredDomains.map(it => `${it}.csv`);
    obj.root.append(ui.divText(`${oneOfStr} the following files should be downloaded to create the view: ${requiredDomains.join(', ')}`));
  }
}