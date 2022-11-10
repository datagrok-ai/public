
import {Observable, Subject} from 'rxjs';
import {IMonomerLib, Monomer} from '../types/index';

export class MonomerLib implements IMonomerLib {
  private _monomers: { [type: string]: { [name: string]: Monomer } } = {};
  private _onChanged = new Subject<any>();

  constructor(monomers: { [type: string]: { [name: string]: Monomer } }) {
    this._monomers = monomers;
  }

  getMonomer(monomerType: string, monomerName: string): Monomer | null {
    if (monomerType in this._monomers! && monomerName in this._monomers![monomerType])
      return this._monomers![monomerType][monomerName];
    else
      return null;
  }

  getTypes(): string[] {
    return Object.keys(this._monomers);
  }

  getMonomerMolsByType(type: string): {[symbol: string]: string} {
    let res: {[symbol: string]: string} = {};

    Object.keys(this._monomers[type]).forEach(monomerSymbol => {
      res[monomerSymbol] = this._monomers[type][monomerSymbol].molfile;
    });

    return res;
  }

  getMonomerNamesByType(type: string): string[] {
    return Object.keys(this._monomers[type]);
  }

  get onChanged(): Observable<any> { 
    return this._onChanged;
  }

  public update(lib: IMonomerLib): void {
    const typesNew = lib.getTypes();
    const types = this.getTypes();

    typesNew.forEach(type => {
      //could possibly rewrite -> TODO: check duplicated monomer symbol

      if (!types.includes(type))
        this._monomers![type] = {};

      const monomers = lib.getMonomerNamesByType(type);
      monomers.forEach(monomerName =>{
        this._monomers[type][monomerName] = lib.getMonomer(type, monomerName)!;
      })
    });

    this._onChanged.next();
  }
}
