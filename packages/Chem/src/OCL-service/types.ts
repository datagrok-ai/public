import * as OCL from 'openchemlib/full';

export interface IChemProperty {
    name: string;
    valueFunc: (mol: OCL.Molecule) => any;
}
