/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {wait} from "datagrok-api/ui";


export let _package = new DG.Package();
let rest = 'https://pubchem.ncbi.nlm.nih.gov/rest';
let pug = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';


//name: PubChem Search Widget
//tags: widgets
//input: string mol {semType: Molecule}
//input: string searchType
//output: widget result
export function PubChemSearchWidget(mol, searchType) {
  let headerHost = ui.divH([]);
  let compsHost = ui.divH([ui.loader()]);
  let panel = ui.divV([compsHost]);
  let search = {
    'similarity': async () => pc.similaritySearch('smiles', mol),
    'substructure': async () => pc.substructureSearch('smiles', mol)
  }

  if (!searchType in search) {
    throw "DrugBankSearch: No such search type" + searchType;
  }

  search[searchType]().then(t => {
      compsHost.removeChild(compsHost.firstChild);
      if (t == null || t.filter.trueCount === 0) {
        compsHost.appendChild(ui.divText('No matches'));
        return;
      }
      t.col('CanonicalSMILES').semType = 'Molecule';
      t.col('CanonicalSMILES').setTag('cell.renderer', 'Molecule');

      compsHost.appendChild(t.plot.grid().root);
      headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
        t.name = `"PubChem Similarity Search"`;
        grok.shell.addTableView(t);
      }, 'Open compounds as table'));
      compsHost.style.overflowY = 'auto';
    }
  ).catch(err => {
    if (compsHost.children.length > 0) {
      compsHost.removeChild(compsHost.firstChild);
    }
    let div = ui.divText('No matches');
    ui.tooltip.bind(div, `${err}`);
    compsHost.appendChild(div);
  });

  return new DG.Widget(panel);
}


//name: PubChem Substructure Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function PubChemSubstructureSearchPanel(mol) {
  return PubChemSearchWidget(mol, 'substructure');
}

//name: PubChem Similarity Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function PubChemSimilaritySearchPanel(mol) {
  return PubChemSearchWidget(mol, 'similarity');
}

class PubChem {
  async similaritySearch(idType, id, params = {}) {
    let listId = await this._asyncSearchId("similarity", idType, id, params)
    let json = undefined;
    let maxRequests = 10;
    while (!json && maxRequests > 0) {
      maxRequests -= 1;
      json = await this._getListById(listId);
    }
    console.error(json);
    let df = DG.DataFrame.fromObjects(json)
    return df

  }

  async identitySearch(idType, id, params = {}) {
    let listId = await this._asyncSearchId("identity", idType, id, params)
    let json = await this._getListById(listId)
    let df = DG.DataFrame.fromObjects(json)
    return df
  }

  async substructureSearch(idType, id, params = {}) {
    let listId = await this._asyncSearchId("substructure", idType, id, params)
    //todo: put in separate function
    let json = undefined;
    let maxRequests = 10;
    while (!json && maxRequests > 0) {
      maxRequests -= 1;
      json = await this._getListById(listId);
    }
    console.error(json);
    let df = DG.DataFrame.fromObjects(json)
    return df
  }

  async smilesToPubChem(smiles) {
    let s = await this.getBy('smiles', 'cids', smiles);
    let cids = s["IdentifierList"]["CID"][0];
    return cids;
  }

  async getBy(idType, idTypeReturn, id, params = {}) {
    return new Promise(function (resolve, reject) {
      let url = `/compound/${idType}/${id}`
      params = params ? params : {}
      let xmlHttp = new XMLHttpRequest();
      let theUrl = `${pug}${url}/${idTypeReturn}/JSON` +
        `${params == {} ? '' : '?' + Object.keys(params).map(key => key + '=' + params[key]).join('&')}`
      xmlHttp.open("GET", theUrl, true);

      xmlHttp.onload = function () {
        if (this.status >= 200 && this.status < 300) {
          let id = JSON.parse(xmlHttp.responseText);
          resolve(id);
        } else {
          reject();
        }
      };
      xmlHttp.onerror = () => reject();
      xmlHttp.send(undefined);
    });
  }

  _getListById(listId, params = {}) {
    return new Promise(function (resolve, reject) {
      let url = `/compound/listkey`
      params = params ? params : {}
      let xmlHttp = new XMLHttpRequest();
      let theUrl = `${pug}${url}/${listId}/property/CanonicalSMILES/JSON` +
        `${params == {} ? '' : '?' + Object.keys(params).map(key => key + '=' + params[key]).join('&')}`
      console.error(theUrl);
      xmlHttp.open("GET", theUrl, true);
      xmlHttp.onload = function () {
        if (this.status >= 200 && this.status < 300) {
          let json = JSON.parse(xmlHttp.responseText);
          console.error(json['PropertyTable']);
          if (json['PropertyTable']) {
            console.error(json['PropertyTable']['Properties']);
            resolve(json['PropertyTable']['Properties']);
          } else {
            resolve(undefined);
          }
        } else
          reject();
      };
      xmlHttp.onerror = () => reject();
      xmlHttp.send(null);
    });
  }

  _asyncSearchId(searchType, idType, id, params = {}) {
    return new Promise(function (resolve, reject) {
      params['MaxRecords'] = 20
      let url = `/compound/${searchType}`
      let xmlHttp = new XMLHttpRequest();
      let theUrl = `${pug}${url}/${idType}/${encodeURIComponent(id)}/JSON` +
        `${params == null ? '' : '?' + Object.keys(params).map(key => key + '=' + params[key]).join('&')}`
      console.error(theUrl)
      xmlHttp.open("GET", theUrl, true);
      xmlHttp.onload = function () {
        if (this.status >= 200 && this.status < 300) {
          let id = JSON.parse(xmlHttp.responseText)['Waiting']['ListKey'];
          resolve(id);
        } else
          reject();
      };
      xmlHttp.onerror = () => reject();
      xmlHttp.send(null);
    });
  }
}

let pc = new PubChem()
