import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable} from 'rxjs';
import {
  HELM_REQUIRED_FIELD as REQ,
  HELM_RGROUP_FIELDS as RGP, HELM_OPTIONAL_FIELDS as OPT, HELM_POLYMER_TYPE
} from '../utils/const';
import {ALPHABET} from '../utils/macromolecule';

export type RGroup = {
  [RGP.CAP_GROUP_SMILES]: string,
  [RGP.ALTERNATE_ID]: string,
  [RGP.CAP_GROUP_NAME]: string,
  [RGP.LABEL]: string,
}
export type Monomer = {
  [REQ.SYMBOL]: string,
  [REQ.NAME]: string,
  [REQ.MOLFILE]: string,
  [REQ.AUTHOR]: string,
  [REQ.ID]: number,
  [REQ.RGROUPS]: RGroup[],
  [REQ.SMILES]: string,
  [REQ.POLYMER_TYPE]: string,
  [REQ.MONOMER_TYPE]: string,
  [REQ.CREATE_DATE]: string | null,
  [OPT.NATURAL_ANALOG]?: string,
  [OPT.META]?: { [property: string]: any },
  lib?: IMonomerLib,
};

export interface IMonomerLib {
  get source(): string | undefined;
  get error(): string | undefined;

  getMonomer(polymerType: string, monomerSymbol: string): Monomer | null;
  getMonomerMolsByPolymerType(polymerType: string): { [monomerSymbol: string]: string } | null;
  getMonomerSymbolsByRGroup(rGroupNumber: number, polymerType: string, element?: string): string[];
  getMonomerSymbolsByType(polymerType: string): string[];
  getPolymerTypes(): string[];
  update(lib: IMonomerLib): void;
  get onChanged(): Observable<any>;

  /** Summary with lib monomer count by type, csv */
  getSummary(): string;
  getTooltip(polymerType: string, monomerSymbol: string): HTMLElement;
}

export const alphabetPolymerTypes = {
  [ALPHABET.DNA]: HELM_POLYMER_TYPE.RNA,
  [ALPHABET.RNA]: HELM_POLYMER_TYPE.RNA,
  [ALPHABET.PT]: HELM_POLYMER_TYPE.PEPTIDE,
  [ALPHABET.UN]: HELM_POLYMER_TYPE.PEPTIDE,
};
