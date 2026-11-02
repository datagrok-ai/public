import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Monomer} from '../types/monomer-library';
import {MonomerType, PolymerType} from '../helm/types';

/** The fields of the monomer editor that a caller reads and sets by name. */
export interface MonomerFormFields {
  /** The structure, as SMILES with its R groups ([*:1], [*:2], ...). */
  molecule: string;
  polymerType: PolymerType;
  monomerType: MonomerType;
  symbol: string;
  name: string;
  id: number | null;
  naturalAnalog: string;
}

/** The monomer editor of the monomer manager. */
export interface INewMonomerForm {
  get form(): HTMLElement;
  /** The fields as the editor holds them now. */
  getFields(): MonomerFormFields;
  /** Sets the fields given; the others keep their values. */
  setFields(fields: Partial<MonomerFormFields>): void;
  setMonomer(monomer: Monomer): void;
}

export interface IMonomerManager {

  /** Creates new monomer library in correct folder and adds given monomers */
  createNewMonomerLib(providerName: string, libName: string, monomers: Monomer[]): Promise<void>;

  editMonomer(dfRow: DG.Row): Promise<Monomer>;

  getNewMonomerForm(): INewMonomerForm;
}
