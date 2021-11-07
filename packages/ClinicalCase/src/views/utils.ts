import { study } from "../clinical-study";
import * as ui from "datagrok-api/ui";
import { domains } from "../sdtm-meta";

export function updateDivInnerHTML(div: HTMLDivElement, content: any){
    div.innerHTML = '';
    div.append(content);
  }

export function checkMissingDomains(requiredDomainsAndCols: any, or: boolean, obj: any) {
  let reqDomains = Object.keys(requiredDomainsAndCols);
  let missingDomains = reqDomains.filter(it => study.domains[it] === null);
  let noMissingDomains = or ? reqDomains.some(it => study.domains[it] !== null) : !missingDomains.length;
  if (noMissingDomains) {
    if (or) {
      reqDomains = reqDomains.filter(it => study.domains[it] !== null);
    }
    if (checkMissingColumns(obj, reqDomains, requiredDomainsAndCols)) {
      obj.createView();
    }
  } else {
    const errorsDiv = ui.div([], { style: { margin: 'auto', textAlign: 'center' } });
    createMissingDataDiv(errorsDiv, missingDomains, 'Missing domains:');
    obj.root.append(errorsDiv);
  }
}

export function checkMissingColumns(obj: any, reqDomains: string[], requiredDomainsAndCols: any) {
  const errorsDiv = ui.divV([], {style: {margin: 'auto', textAlign: 'center'}});
  let noMissingCols = true;
  reqDomains.forEach(domain => {
    const domainColumns = study.domains[domain].columns.names();
    const missingColumns = requiredDomainsAndCols[domain].filter(it => !domainColumns.includes(it));
    if(missingColumns.length){
      noMissingCols = false;
      createMissingDataDiv(errorsDiv, missingColumns, `Missing columns in ${domain}:`)
    }
  })
  if(!noMissingCols){
    obj.root.append(errorsDiv);
  }
  return noMissingCols;
}

function createMissingDataDiv(div: HTMLDivElement, missingDomainsOrCols: string[], header: string){
  let domainsDiv = ui.div();
  missingDomainsOrCols.forEach(it => {domainsDiv.append(ui.divText(it))})
  div.append(ui.div([
    ui.h2(header),
    domainsDiv
  ]));
}