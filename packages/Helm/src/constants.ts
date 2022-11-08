export const jsonSdfMonomerLibDict = {
    "monomerType": null, 
    "smiles": null, 
    "name": "MonomerName", 
    "author": null, 
    "molfile": "molecule", 
    "naturalAnalog": "MonomerNaturalAnalogCode", 
    "rgroups": "MonomerCaps",
    "createDate": null, 
    "id": null, 
    "polymerType": "MonomerType", 
    "symbol": "MonomerCode"
}

export type WebEditorMonomer = {
    id: string,     //symbol
    n: string,      //name
    na: string,     //natural analog
    type: string,   //polymer type
    mt: string,     //monomer type
    m: string,      //molfile
    at: {[group: string]: string},   //substituents
    rs: number      //number of substituents
  };

export const SMILES = 'smiles';
export const RGROUPS = "rgroups";
export const MONOMER_SYMBOL = "symbol";
export const RGROUP_CAP_GROUP_SMILES = "capGroupSmiles";
export const RGROUP_ALTER_ID = "alternateId";
export const RGROUP_CAP_GROUP_NAME = "capGroupName";
export const RGROUP_LABEL = "label";
export const SDF_MONOMER_NAME = "MonomerName";