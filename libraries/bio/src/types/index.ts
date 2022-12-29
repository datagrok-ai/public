import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable} from 'rxjs';


export type Monomer = {
  symbol: string,
  name: string,
  naturalAnalog: string,
  molfile: string,
  polymerType: string,
  monomerType: string,
  rgroups: {capGroupSmiles: string, alternateId: string, capGroupName: string, label: string }[],
  data: {[property: string]: any}
};

export interface IMonomerLib {
  getMonomer(monomerType: string, monomerName: string): Monomer | null;
  getMonomerMolsByType(type: string): {[symbol: string]: string} | null;
  getMonomerNamesByType(type: string): string[];
  getTypes(): string[];
  update(lib: IMonomerLib): void;
  get onChanged(): Observable<any>;
}
