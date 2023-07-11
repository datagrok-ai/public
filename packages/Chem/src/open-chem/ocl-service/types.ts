import * as OCL from 'openchemlib/full';

export interface IChemProperty {
    name: string;
    valueFunc: (mol: OCL.Molecule) => any;
    type: IChemPropertyType;
}
// we need this instead of DG.TYPE.FLOAT or int, because anything used by workers can not import DG, as it is external
export type IChemPropertyType = 'float' | 'int';
