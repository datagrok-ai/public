import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Monomer} from '../types/monomer-library';
import {MonomerTypes, PolymerTypes} from '../helm/consts';
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

export interface IMonomerGallery {
  get monomerGallery(): HTMLElement;
  groupBy(by: string): void;
  filterBySearch(search: string): void;
}

export interface IMonomerManager {

  /** Creates new monomer library in correct folder and adds given monomers */
  createNewMonomerLib(providerName: string, libName: string, monomers: Monomer[]): Promise<void>;

  editMonomer(dfRow: DG.Row): Promise<Monomer>;

  getNewMonomerForm(): INewMonomerForm;

  //getMonomerGallery(): IMonomerGallery;
}

export const MonomerInputProperties: { [key: string]: DG.Property } = {
  'monomerType': DG.Property.js('monomerType', DG.TYPE.STRING,
    {choices: [MonomerTypes.BACKBONE, MonomerTypes.BRANCH, MonomerTypes.TERMINAL]}),
  'molecule': DG.Property.fromOptions({name: 'molecule', type: DG.InputType.Molecule}),
  'name': DG.Property.js('name', DG.TYPE.STRING),
  'naturalAnalog': DG.Property.js('naturalAnalog', DG.TYPE.STRING),
  'id': DG.Property.js('id', DG.TYPE.INT),
  'polymerType': DG.Property.js('polymerType', DG.TYPE.STRING,
    {choices: [PolymerTypes.RNA, PolymerTypes.PEPTIDE, PolymerTypes.CHEM, PolymerTypes.BLOB, PolymerTypes.G]}),
  'symbol': DG.Property.js('symbol', DG.TYPE.STRING),
};
