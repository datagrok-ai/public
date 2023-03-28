import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable} from 'rxjs';
import {HELM_MONOMER_TYPE, HELM_POLYMER_TYPE} from '../utils/const';

export type Monomer = {
  symbol: string,
  name: string,
  naturalAnalog: string,
  molfile: string,
  polymerType: HELM_POLYMER_TYPE,
  monomerType: HELM_MONOMER_TYPE,
  rgroups: {capGroupSmiles: string, alternateId: string, capGroupName: string, label: string }[],
  data: {[property: string]: any}
};

export interface IMonomerLib {
  getMonomer(polymerType: HELM_POLYMER_TYPE, monomerName: string): Monomer | null;
  getMonomerMolsByType(polymerType: HELM_POLYMER_TYPE): {[symbol: string]: string} | null;
  getMonomerNamesByType(polymerType: HELM_POLYMER_TYPE): string[];
  getTypes(): string[];
  update(lib: IMonomerLib): void;
  get onChanged(): Observable<any>;
}
