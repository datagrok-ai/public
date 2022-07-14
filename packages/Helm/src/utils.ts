import * as DG from 'datagrok-api/dg';
import { RGROUP_CAP_GROUP_NAME, RGROUP_CAP_GROUP_SMILES, jsonSdfMonomerLibDict, MONOMER_SYMBOL, RGROUP_ALTER_ID, RGROUPS, RGROUP_LABEL, SDF_MONOMER_NAME } from "./constants";

export function createJsonMonomerLibFromSdf(table: DG.DataFrame): any {
    const resultLib = [];
    for (let i = 0; i < table.rowCount; i++) {
      const monomer: { [key: string]: string | any } = {};
      Object.keys(jsonSdfMonomerLibDict).forEach(key => {
        if (key === MONOMER_SYMBOL) {
          const monomerSymbol = table.get(jsonSdfMonomerLibDict[key], i);
          monomer[key] = monomerSymbol === '.' ? table.get(SDF_MONOMER_NAME, i) : monomerSymbol;
        } else if (key === RGROUPS) {
          const rgroups = table.get(jsonSdfMonomerLibDict[key], i).split('\n');
          const jsonRgroups: any[] = [];
          rgroups.forEach((g: string) => {
            const rgroup: { [key: string]: string | any } = {};
            const altAtom = g.substring(g.lastIndexOf("]") + 1);
            let radicalNum = g.match(/\[R(\d+)\]/)![1];
            rgroup[RGROUP_CAP_GROUP_SMILES] = altAtom === 'H' ? `[*:${radicalNum}][H]` : `O[*:${radicalNum}]`;
            rgroup[RGROUP_ALTER_ID] = altAtom === 'H' ? `R${radicalNum}-H` : `R${radicalNum}-OH`;
            rgroup[RGROUP_CAP_GROUP_NAME] = altAtom === 'H' ? `H` : `OH`;
            rgroup[RGROUP_LABEL] = `R${radicalNum}`;
            jsonRgroups.push(rgroup);
          })
          monomer[key] = jsonRgroups;
        } else {
          if((jsonSdfMonomerLibDict as { [key: string]: string | any })[key]) {
            monomer[key] = table.get((jsonSdfMonomerLibDict as { [key: string]: string | any })[key], i);
          }
        }
      })
      resultLib.push(monomer);
    }
    return resultLib;
  }