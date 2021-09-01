/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {JOIN_TYPE} from "datagrok-api/dg";

export let _package = new DG.Package();


//name: DrugBankSubstructureSearch
//input: string mol {semType: Molecule}
//input: bool searchType
//output: dataframe dbResult
export async function
drugBankSubstructureSearch(mol, substructLibrary) {


    try {
        let dbdf = await grok.data.loadTable(`${_package.webRoot}db.csv`);

        if (dbdf == null) {
            return DG.DataFrame.create();
        }
        let bitset = await grok.chem.substructureSearch(dbdf.getCol('molecule'), mol, {'substructLibrary': substructLibrary})
        if (bitset != null) {
            dbdf.filter.copyFrom(bitset);
            return dbdf;
        }

        return null;
    } catch (e) {
        console.error("In similarityScoring: " + e.toString());
        throw e;
    }
}

//name: DrugBankSimilaritySearch
//input: string molecule {semType: Molecule}
//input: int limit
//input: int cutoff
//output: dataframe drbResult
export async function drugBankSimilaritySearch(molecule, limit, cutoff) {
    try {
        let dbdf = await grok.data.loadTable(`${_package.webRoot}db.csv`);

        if (dbdf == null) {
            return DG.DataFrame.create();
        }
        let searchdf = await grok.chem.findSimilar(
            dbdf.getCol('molecule'), molecule, {'limit': limit, 'cutoff': cutoff});
        if (searchdf == null) {
            return DG.DataFrame.create();
        }
        let index = DG.Column.fromInt32Array('index', Int32Array.from(new Array(dbdf.rowCount), (x, i) => i))
        dbdf.columns.add(index)
        return dbdf.join(searchdf, ['index'], ['index'], ['molecule'], ['molecule'], JOIN_TYPE.INNER, true)


    } catch (e) {
        console.error(e);
        return null;
    }

}

//name: DrugBank Search Widget
//tags: widgets
//input: string mol {semType: Molecule}
//input: string searchType
//output: widget result
export function DrugBankSearchWidget(mol, searchType) {
    let headerHost = ui.divH([]);
    let compsHost = ui.divH([ui.loader()]);
    let panel = ui.divV([compsHost]);
    let search = {
        'similarity': async () => drugBankSimilaritySearch(mol, 20, 0),
        'substructure': async () => drugBankSubstructureSearch(mol, false)
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

            function getTooltip(n) {
                let props = {
                    'DRUGBANK_ID': t.get('DRUGBANK_ID', n),
                };
                return ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]);
            }

            let piv = -1;
            import('openchemlib/full.js').then((OCL) => {
                for (let n = 0; n < Math.min(t.rowCount, 20); n++) {
                    piv = t.filter.findNext(piv + 1, true)


                    if (piv < 0) {
                        break;
                    }

                    let smiles = t.get('molecule', piv);

                    let molecule = document.createElement('div');

                    let m = OCL.Molecule.fromMolfile(smiles);
                    molecule.innerHTML = m.toSVG(250, 150);

                    ui.tooltip.bind(molecule, () => getTooltip(n));
                    molecule.addEventListener('click', function () {
                        window.open(t.get('link', piv), '_blank');
                    });
                    compsHost.appendChild(molecule);
                }
            });
            headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
                t.name = `"DrugBank Similarity Search"`;
                grok.shell.addTableView(t);
            }, 'Open compounds as table'));
            compsHost.style.overflowY = 'auto';
        }
    )
        .catch(err => {
            if (compsHost.children.length > 0) {
                compsHost.removeChild(compsHost.firstChild);
            }
            let div = ui.divText('No matches');
            ui.tooltip.bind(div, `${err}`);
            compsHost.appendChild(div);
        });
    return new DG.Widget(panel);
}


//name: DrugBank Substructure Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function DrugBankSubstructureSearchPanel(mol) {
    return DrugBankSearchWidget(mol, 'substructure');
}

//name: DrugBank Similarity Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function DrugBankSimilaritySearchPanel(mol) {

    return DrugBankSearchWidget(mol, 'similarity');
}




