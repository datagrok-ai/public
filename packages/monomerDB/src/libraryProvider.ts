/* eslint-disable max-len */
import {IMonomerLibProvider, IMonomerSet, Monomer, MonomerLibData} from '@datagrok-libraries/bio/src/types/monomer-library';
import * as api from './package-api';
import {PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {Observable, Subject} from 'rxjs';
export class DBLibraryProvider implements IMonomerLibProvider {
  name = 'DBLibrary';
  private isLoadedOnce: boolean = false;
  private listsOfLibraries: string[] = [];
  private _loadPromise: Promise<MonomerLibData[]> = Promise.resolve([]);

  private _onChanged: Subject<void> = new Subject<void>();
  get onChanged(): Observable<void> {
    return this._onChanged;
  }

  get loadPromise() {
    return this._loadPromise.then(() => { });
  }
  async listLibraries(): Promise<string[]> {
    if (!this.isLoadedOnce)
      await this.refreshLists();
    return this.listsOfLibraries;
  }
  async refreshLists(): Promise<void> {
    const newList = (await api.queries.listLibraries())?.col('name')?.toList() ?? [];

    if (newList.length !== this.listsOfLibraries.length ||
      !newList.every((value, index) => value === this.listsOfLibraries[index])) {
      this.listsOfLibraries = newList;
      this._onChanged.next();
    }
    this.isLoadedOnce = true;
  }

  async loadLibraries(names: string[]): Promise<MonomerLibData[]> {
    this._loadPromise = this._loadPromise.then(async () => {
      const monomerData: MonomerLibData[] = [];
      for (const name of names) {
        const df = await api.queries.getMonomersByLibrary(name);
        if (!df) {
          console.warn(`Library ${name} not found in DB`);
          continue;
        }
        const data: MonomerLibData = {};
        for (const row of df.rows) {
          const monomer: Monomer = {
            id: row.get('id'),
            symbol: row.get('symbol'),
            name: row.get('name'),
            monomerType: row.get('monomerType'),
            polymerType: row.get('polymerType'),
            smiles: row.get('smiles'),
            molfile: row.get('molfile'),
            rgroups: JSON.parse(row.get('rgroups')),
            naturalAnalog: row.get('naturalAnalog'),
            author: row.get('author'),
            createDate: row.get('createDate'),
            meta: JSON.parse(row.get('meta') || '{}'),
          };
          if (!data[monomer.polymerType])
            data[monomer.polymerType] = {};
          data[monomer.polymerType][monomer.symbol] = monomer;
        }
        monomerData.push(data);
      }


      return monomerData;
    });
    return this._loadPromise;
  }

  async deleteLibrary(name: string): Promise<void> {
    //not implemented
    return Promise.resolve();
  }
  async addOrUpdateLibrary(libraryName: string, monomers: Monomer[]): Promise<void> {
    //not implemented
    return Promise.resolve();
  }

  async addOrUpdateLibraryString(name: string, contentString: string): Promise<void> {
    await this.addOrUpdateLibrary(name, JSON.parse(contentString) as Monomer[]);
    return Promise.resolve();
  }

  getRgroupNumbers(monomer: Monomer): (number | null)[] {
    // TODO: implement proper parsing of R-groups from DB
    const rgroupNumbers: (number | null)[] = new Array(8).fill(null);
    let counter = 0;
    for (const rgroup of monomer.rgroups) {
      if (rgroup.capGroupName?.toLowerCase() === 'h')
        rgroupNumbers[counter++] = 1;
      else
        rgroupNumbers[counter++] = 2;
    }
    return rgroupNumbers;
  }

  async updateOrAddMonomersInLibrary(libraryName: string, monomers: Monomer[]): Promise<void> {
    for (const monomer of monomers) {
      const rgroupNumbers = this.getRgroupNumbers(monomer);
      await api.queries.upsertMonomer(monomer.symbol, monomer.name, monomer.monomerType, monomer.smiles, monomer.molfile,
        libraryName, 'Admin', monomer.polymerType, '', monomer.naturalAnalog ?? null, rgroupNumbers[0], rgroupNumbers[1], rgroupNumbers[2],
        rgroupNumbers[3], rgroupNumbers[4], rgroupNumbers[5], rgroupNumbers[6], rgroupNumbers[7]);
    }
  }

  getSingleLibrary(name: string): Promise<MonomerLibData | null> {
    return this.loadLibraries([name]).then((data) => data.length > 0 ? data[0] : null);
  }

  getLibraryAsString(libName: string): Promise<string> {
    return this.getSingleLibrary(libName).then((data) => data ? JSON.stringify(data) : '');
  }

  listSets(): Promise<string[]> {
    return Promise.resolve([]); // Not implemented
  }

  loadSets(names: string[]): Promise<IMonomerSet[]> {
    return Promise.resolve([]); // Not implemented
  }
  deleteSet(name: string): Promise<void> {
    return Promise.resolve(); // Not implemented
  }

  addOrUpdateSet(setName: string, monomerSet: IMonomerSet): Promise<void> {
    return Promise.resolve(); // Not implemented
  }
  addOrUpdateSetString(name: string, contentString: string): Promise<void> {
    return Promise.resolve(); // Not implemented
  }

  deleteMonomersFromLibrary(libraryName: string, monomers: ({ polymerType: PolymerType; symbol: string; }[])): Promise<void> {
    return Promise.resolve(); // Not implemented
  }
}
