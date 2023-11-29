/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable, Subject} from 'rxjs';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/index';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';

import {_package} from '../../package';

/** Wrapper for monomers obtained from different sources. For managing monomere
 * libraries, use MolfileHandler class instead */
export class MonomerLib implements IMonomerLib {
  public readonly error: string | undefined;

  private _monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } } = {};
  private _onChanged = new Subject<any>();

  constructor(monomers: { [polymerType: string]: { [monomerSymbol: string]: Monomer } }, error?: string) {
    this._monomers = monomers;
    this.error = error;
  }

  getMonomer(polymerType: string, monomerSymbol: string): Monomer | null {
    if (polymerType in this._monomers! && monomerSymbol in this._monomers![polymerType])
      return this._monomers![polymerType][monomerSymbol];
    else
      return null;
  }

  getPolymerTypes(): string[] {
    return Object.keys(this._monomers);
  }

  getMonomerMolsByPolymerType(polymerType: string): { [monomerSymbol: string]: string } {
    const res: { [monomerSymbol: string]: string } = {};

    Object.keys(this._monomers[polymerType] ?? {}).forEach((monomerSymbol) => {
      res[monomerSymbol] = this._monomers[polymerType][monomerSymbol].molfile;
    });

    return res;
  }

  getMonomerSymbolsByType(polymerType: string): string[] {
    return Object.keys(this._monomers[polymerType]);
  }

  /** Get a list of monomers with specified element attached to specified
   * R-group
   * WARNING: RGroup numbering starts from 1, not 0*/
  getMonomerSymbolsByRGroup(rGroupNumber: number, polymerType: string, element?: string): string[] {
    const monomerSymbols = this.getMonomerSymbolsByType(polymerType);
    let monomers = monomerSymbols.map((sym) => this.getMonomer(polymerType, sym));
    monomers = monomers.filter((el) => el !== null);
    if (monomers.length === 0)
      return [];

    function findAllIndices<T>(arr: T[], element: T): number[] {
      return arr.map((value, index) => (value === element ? index : -1))
        .filter((index) => index !== -1);
    }

    monomers = monomers.filter((monomer) => {
      if (!monomer?.rgroups)
        return false;
      let criterion = monomer?.rgroups.length >= rGroupNumber;
      const molfileHandler = MolfileHandler.getInstance(monomer.molfile);
      const rGroupIndices = findAllIndices(molfileHandler.atomTypes, 'R#');
      criterion &&= true;
      return criterion;
    });
    return monomers.map((monomer) => monomer?.symbol!);
  }

  get onChanged(): Observable<any> {
    return this._onChanged;
  }

  private _updateInt(lib: IMonomerLib): void {
    const typesNew = lib.getPolymerTypes();
    const types = this.getPolymerTypes();

    typesNew.forEach((type) => {
      //could possibly rewrite -> TODO: check duplicated monomer symbol

      if (!types.includes(type))
        this._monomers![type] = {};

      const monomers = lib.getMonomerSymbolsByType(type);
      monomers.forEach((monomerSymbol) => {
        this._monomers[type][monomerSymbol] = lib.getMonomer(type, monomerSymbol)!;
      });
    });
  }

  public update(lib: IMonomerLib): void {
    this._updateInt(lib);
    this._onChanged.next();
  }

  public updateLibs(libList: IMonomerLib[], reload: boolean = false): void {
    if (reload) this._monomers = {};
    for (const lib of libList) if (!lib.error) this._updateInt(lib);
    this._onChanged.next();
  }

  public clear(): void {
    this._monomers = {};
    this._onChanged.next();
  }
}
