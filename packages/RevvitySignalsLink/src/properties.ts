import * as DG from 'datagrok-api/dg';

export const created = DG.Property.create('Created', DG.TYPE.DATE_TIME, (x: any) => x, (x: any, v) => x = v);
export const edited = DG.Property.create('Edited', DG.TYPE.DATE_TIME, (x: any) => x, (x: any, v) => x = v);
export const creator = DG.Property.create('Creator', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
export const editor = DG.Property.create('Editor', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
export const structure = DG.Property.create('Structure', DG.TYPE.STRING, (x: any) => x, (x: any, v) => x = v);
export const id = DG.Property.create('Id_float', DG.TYPE.FLOAT, (x: any) => x, (x: any, v) => x = v);
export const idInt = DG.Property.create('Id_int', DG.TYPE.INT, (x: any) => x, (x: any, v) => x = v);
structure.semType = DG.SEMTYPE.MOLECULE;

export function getDefaultProperties(): DG.Property[] {
    return [created, edited, creator, editor, structure, id, idInt];
}


const names = ['Davit Rizhinashvili', 'Maria Dolotova', 'Andrew Skalkin', 'Ed Jaeger'];

function createdBySuggestions(text: string) {
    return names.filter((it) => it.toLowerCase().includes(text.toLowerCase()));
}