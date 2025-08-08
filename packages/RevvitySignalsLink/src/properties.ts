import * as DG from 'datagrok-api/dg';

export const created = DG.Property.create('Created', DG.TYPE.DATE_TIME, (x: any) => x, (x: any, v) => x = v);
export const createdBy = DG.Property.create('Created by', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
createdBy.options.suggestionsFunction = createdBySuggestions;
export const id = DG.Property.create('Id', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
export const mw = DG.Property.create('MW', DG.TYPE.FLOAT, (x: any) => x, (x: any, v) => x = v);
export const molecule = DG.Property.create('Molecule', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
molecule.semType = DG.SEMTYPE.MOLECULE;
// export const chem = DG.Property.create('Chem', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
// chem.semType = DG.SEMTYPE.MOLECULE;

export function getProperties(): DG.Property[] {
    return [created, createdBy, id, mw, molecule];
}


const names = ['Davit Rizhinashvili', 'Maria Dolotova', 'Andrew Skalkin', 'Ed Jaeger'];

function createdBySuggestions(text: string) {
    return names.filter((it) => it.toLowerCase().includes(text.toLowerCase()));
}