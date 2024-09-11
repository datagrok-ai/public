import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Monomer} from '../types/index';
import {MonomerTypes, PolymerTypes} from '../helm/consts';

export interface INewMonomerForm {
  get fieldInputs(): { [key: string]: DG.InputBase<any> | grok.chem.Sketcher };
  get form(): HTMLElement;
  get rgroupInputs(): { [key: string]: DG.InputBase<any> }[];
  get metaInputs(): { [key: string]: DG.InputBase<any> }[];
  setMonomer(monomer: Monomer): void;
}

export interface IMonomerGallery {
  get monomerGallery(): HTMLElement;
  groupBy(by: string): void;
  filterBySearch(search: string): void;
}

export interface IMonomerManager {

  /** Creates new monomer library in correct folder and adds given monomers */
  createNewMonomerLib(libName: string, monomers: Monomer[]): Promise<void>;

  cloneMonomer(dfRow: DG.Row): Monomer;

  getNewMonomerForm(): INewMonomerForm;

  //getMonomerGallery(): IMonomerGallery;
}

export const MonomerInputProperties: { [key: string]: DG.Property } = {
  'monomerType': DG.Property.js('monomerType', DG.TYPE.STRING, {choices: [MonomerTypes.BACKBONE, MonomerTypes.BRANCH, MonomerTypes.TERMINAL]}),
  'molecule': DG.Property.fromOptions({name: 'molecule', type: DG.InputType.Molecule}),
  'name': DG.Property.js('name', DG.TYPE.STRING),
  'naturalAnalog': DG.Property.js('naturalAnalog', DG.TYPE.STRING),
  'id': DG.Property.js('id', DG.TYPE.INT),
  'polymerType': DG.Property.js('polymerType', DG.TYPE.STRING, {choices: [PolymerTypes.RNA, PolymerTypes.PEPTIDE, PolymerTypes.CHEM, PolymerTypes.BLOB, PolymerTypes.G]}),
  'symbol': DG.Property.js('symbol', DG.TYPE.STRING),
};
