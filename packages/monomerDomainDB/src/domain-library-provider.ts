/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';

import {
  IMonomerLib, IMonomerLibProvider, IMonomerLinkData, IMonomerSet, Monomer, MonomerLibData,
  MonomerSet, MonomerSetPlaceholder, getMonomerLibHelper,
} from '@datagrok-libraries/bio/src/types/monomer-library';
import {MonomerType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {Observable, Subject} from 'rxjs';

import {
  MonomerColumn, MonomerInsert, MonomerRow, SetPlaceholderInsert, SetPlaceholderRow, monomersDb,
} from './generated/db';

export const DOMAIN_DB_PROVIDER_NAME = 'DomainDB';

const PAGE_SIZE = 1000;

interface DeletableClient {
  deleteWhere(filter: any, options?: {limit?: number}): Promise<DG.DomainDeleteReport>;
}

/** Monomer library provider backed by the `monomers` domain schema: libraries and monomers are
 * entity-mapped rows with server-enforced security, audit, optimistic concurrency, and
 * business-key upsert — instead of whole-file JSON rewrites. */
export class DomainDbLibraryProvider implements IMonomerLibProvider {
  name: string = DOMAIN_DB_PROVIDER_NAME;

  private static _instance: DomainDbLibraryProvider | null = null;
  static getInstance(): DomainDbLibraryProvider {
    if (DomainDbLibraryProvider._instance == null) {
      const inst = DomainDbLibraryProvider._instance = new DomainDbLibraryProvider();
      inst.enqueue(() => inst.refreshLists());
    }
    return DomainDbLibraryProvider._instance;
  }

  private libNames: string[] = [];
  private setNames: string[] = [];
  private listsLoaded: boolean = false;

  private _queue: Promise<void> = Promise.resolve();
  private readonly _onChanged = new Subject<void>();

  get onChanged(): Observable<void> { return this._onChanged; }
  get loadPromise(): Promise<void> { return this._queue; }

  /** Serializes mutating operations so overlapping calls cannot interleave. */
  private enqueue<T>(fn: () => Promise<T>): Promise<T> {
    const run = this._queue.then(fn);
    this._queue = run.then(() => undefined, () => undefined);
    return run;
  }

  // -- Libraries --

  async listLibraries(): Promise<string[]> {
    if (!this.listsLoaded)
      await this.refreshLists();
    return this.libNames;
  }

  refreshLists(): Promise<void> {
    return this.reloadLists(false);
  }

  async loadLibraries(names: string[]): Promise<MonomerLibData[]> {
    const libs: MonomerLibData[] = [];
    for (const name of names.map(normName)) {
      let data: MonomerLibData = {};
      try {
        const libRow = await this.findLibrary(name);
        if (libRow == null) {
          console.warn(`DomainDbLibraryProvider: library '${name}' not found`);
          continue;
        }
        data = await this.buildLibData(libRow.id);
      } catch (err) {
        console.error(`DomainDbLibraryProvider: failed to load library '${name}'`, err);
        continue;
      } finally {
        libs.push(data);
      }
    }
    return libs;
  }

  async getSingleLibrary(name: string): Promise<MonomerLibData | null> {
    const libRow = await this.findLibrary(normName(name));
    return libRow == null ? null : await this.buildLibData(libRow.id);
  }

  async getLibraryAsString(libName: string): Promise<string> {
    const name = normName(libName);
    const libRow = await this.findLibrary(name);
    if (libRow == null)
      throw new Error(`Monomer library '${name}' not found`);
    const monomers = (await this.queryLibraryMonomers(libRow.id)).map(rowToMonomer);
    return JSON.stringify(monomers, null, 2);
  }

  async addOrUpdateLibraryString(name: string, contentString: string): Promise<void> {
    const raw = JSON.parse(contentString);
    if (!Array.isArray(raw))
      throw new Error(`Monomer library '${name}' content must be a JSON array of monomers`);
    for (let i = 0; i < raw.length; i++) {
      if (!raw[i]?.symbol || !raw[i]?.polymerType)
        throw new Error(`Monomer library '${name}': entry ${i} lacks 'symbol' or 'polymerType'`);
    }
    await this.addOrUpdateLibrary(name, raw as Monomer[]);
  }

  /** Full-set replace: monomers absent from [monomers] are removed from the library
   * (same semantics as overwriting a library file). */
  addOrUpdateLibrary(libraryName: string, monomers: Monomer[]): Promise<void> {
    const name = normName(libraryName);
    return this.enqueue(async () => {
      const {id: libId} = await monomersDb.libraries.upsert({name});
      await this.upsertMonomers(libId, monomers);
      const keep = new Set(monomers.map((m) => monomerKey(m.polymerType, m.symbol)));
      const stale = (await this.queryLibraryMonomers(libId, ['symbol', 'polymer_type']))
        .filter((r) => !keep.has(monomerKey(r.polymer_type, r.symbol)))
        .map((r) => r.id);
      await this.deleteByIds(monomersDb.monomers, stale);
      await this.reloadLists(true);
    });
  }

  /** Merge: existing monomers (by polymer type + symbol) are updated, new ones appended. */
  updateOrAddMonomersInLibrary(libraryName: string, monomers: Monomer[]): Promise<void> {
    const name = normName(libraryName);
    return this.enqueue(async () => {
      const libRow = await this.findLibrary(name);
      if (libRow == null)
        throw new Error(`Monomer library '${name}' not found`);
      await this.upsertMonomers(libRow.id, monomers);
      await this.reloadLists(true);
    });
  }

  deleteMonomersFromLibrary(libraryName: string,
    monomers: ({polymerType: PolymerType, symbol: string}[])): Promise<void> {
    const name = normName(libraryName);
    return this.enqueue(async () => {
      const libRow = await this.findLibrary(name);
      if (libRow == null)
        return;
      const toDelete = new Set(monomers.map((m) => monomerKey(m.polymerType, m.symbol)));
      const ids = (await this.queryLibraryMonomers(libRow.id, ['symbol', 'polymer_type']))
        .filter((r) => toDelete.has(monomerKey(r.polymer_type, r.symbol)))
        .map((r) => r.id);
      await this.deleteByIds(monomersDb.monomers, ids);
      await this.reloadLists(true);
    });
  }

  deleteLibrary(name: string): Promise<void> {
    const libName = normName(name);
    return this.enqueue(async () => {
      const libRow = await this.findLibrary(libName);
      if (libRow == null)
        return;
      await monomersDb.libraries.delete(libRow.id);
      await this.reloadLists(true);
    });
  }

  // -- Sets --

  async listSets(): Promise<string[]> {
    if (!this.listsLoaded)
      await this.refreshLists();
    return this.setNames;
  }

  async loadSets(names: string[]): Promise<IMonomerSet[]> {
    const sets: IMonomerSet[] = [];
    const lib = await this.tryGetMonomerLib();
    for (const name of names.map(normName)) {
      try {
        const setRow = await monomersDb.monomerSets.first({filter: DG.cond('name', '=', name)});
        if (setRow == null)
          continue;
        const rows = await this.queryAll<SetPlaceholderRow>((limit, offset) =>
          monomersDb.setPlaceholders.query({filter: DG.cond('set_id', '=', setRow.id), sort: 'symbol', limit, offset}));
        const placeholders = rows.map((r) => new MonomerSetPlaceholder(
          lib, r.symbol, r.polymer_type as PolymerType, (r.monomer_type ?? 'Undefined') as MonomerType,
          parseJson<IMonomerLinkData[]>(r.monomer_links, [])));
        sets.push(new MonomerSet(setRow.description ?? '', placeholders, name));
      } catch (err) {
        console.error(`DomainDbLibraryProvider: failed to load monomer set '${name}'`, err);
        continue;
      }
    }
    return sets;
  }

  addOrUpdateSet(setName: string, monomerSet: IMonomerSet): Promise<void> {
    return this.writeSet(normName(setName), monomerSet.description,
      monomerSet.placeholders.map((ph) => ({
        symbol: ph.symbol,
        polymer_type: ph.polymerType as SetPlaceholderInsert['polymer_type'],
        monomer_type: ph.monomerType as SetPlaceholderInsert['monomer_type'],
        monomer_links: JSON.stringify(ph.monomerLinks),
      })));
  }

  /** Accepts the monomer-set file format:
   * `{description, placeholders: {SYMBOL: {polymerType, monomerType, set: [links]}}}`. */
  addOrUpdateSetString(name: string, contentString: string): Promise<void> {
    const raw = JSON.parse(contentString);
    const placeholders = Object.entries(raw?.placeholders ?? {}).map(([symbol, v]: [string, any]) => ({
      symbol,
      polymer_type: v?.polymerType as SetPlaceholderInsert['polymer_type'],
      monomer_type: v?.monomerType as SetPlaceholderInsert['monomer_type'],
      monomer_links: JSON.stringify(v?.set ?? []),
    }));
    return this.writeSet(normName(name), raw?.description ?? '', placeholders);
  }

  deleteSet(name: string): Promise<void> {
    const setName = normName(name);
    return this.enqueue(async () => {
      const setRow = await monomersDb.monomerSets.first({filter: DG.cond('name', '=', setName)});
      if (setRow == null)
        return;
      await monomersDb.monomerSets.delete(setRow.id);
      await this.reloadLists(true);
    });
  }

  // -- Internals --

  private findLibrary(name: string) {
    return monomersDb.libraries.first({filter: DG.cond('name', '=', name)});
  }

  // system columns cannot be NAMED in `columns`, but are always returned alongside the projection
  private queryLibraryMonomers(libraryId: string, columns?: MonomerColumn[]): Promise<MonomerRow[]> {
    return this.queryAll<MonomerRow>((limit, offset) => monomersDb.monomers.query({
      filter: DG.cond('library_id', '=', libraryId), sort: 'symbol',
      ...(columns == null ? {} : {columns}), limit, offset,
    }));
  }

  private async buildLibData(libraryId: string): Promise<MonomerLibData> {
    const data: MonomerLibData = {};
    for (const row of await this.queryLibraryMonomers(libraryId)) {
      const monomer = rowToMonomer(row);
      (data[monomer.polymerType] ??= {})[monomer.symbol] = monomer;
    }
    return data;
  }

  private async upsertMonomers(libraryId: string, monomers: Monomer[]): Promise<void> {
    if (monomers.length === 0)
      return;
    const report = await monomersDb.monomers.batch(monomersToDf(libraryId, monomers), {mode: 'upsert'});
    checkBatchReport(report, 'monomers');
  }

  private writeSet(name: string, description: string,
    placeholders: Omit<SetPlaceholderInsert, 'set_id'>[]): Promise<void> {
    return this.enqueue(async () => {
      const {id: setId} = await monomersDb.monomerSets.upsert({name, description});
      if (placeholders.length > 0) {
        const report = await monomersDb.setPlaceholders.batch(
          placeholders.map((ph) => ({...ph, set_id: setId})), {mode: 'upsert'});
        checkBatchReport(report, 'set placeholders');
      }
      const keep = new Set(placeholders.map((ph) => ph.symbol));
      const stale = (await this.queryAll<SetPlaceholderRow>((limit, offset) =>
        monomersDb.setPlaceholders.query({filter: DG.cond('set_id', '=', setId), limit, offset})))
        .filter((r) => !keep.has(r.symbol)).map((r) => r.id);
      await this.deleteByIds(monomersDb.setPlaceholders, stale);
      await this.reloadLists(true);
    });
  }

  private async deleteByIds(client: DeletableClient, ids: string[]): Promise<void> {
    if (ids.length === 0)
      return;
    for (;;) {
      const report = await client.deleteWhere(DG.cond('id', '=', ids));
      if (!report.hasMore)
        return;
    }
  }

  private async queryAll<T>(page: (limit: number, offset: number) => Promise<T[]>): Promise<T[]> {
    const res: T[] = [];
    for (let offset = 0; ; offset += PAGE_SIZE) {
      const rows = await page(PAGE_SIZE, offset);
      res.push(...rows);
      if (rows.length < PAGE_SIZE)
        return res;
    }
  }

  private async queryNames(page: (limit: number, offset: number) => Promise<{name: string}[]>): Promise<string[]> {
    return (await this.queryAll(page)).map((r) => r.name);
  }

  private async reloadLists(fireAlways: boolean): Promise<void> {
    const [libNames, setNames] = await Promise.all([
      this.queryNames((limit, offset) =>
        monomersDb.libraries.query({columns: ['name'], sort: 'name', limit, offset})),
      this.queryNames((limit, offset) =>
        monomersDb.monomerSets.query({columns: ['name'], sort: 'name', limit, offset})),
    ]);
    const changed = !sameList(libNames, this.libNames) || !sameList(setNames, this.setNames);
    this.libNames = libNames;
    this.setNames = setNames;
    this.listsLoaded = true;
    if (fireAlways || changed)
      this._onChanged.next();
  }

  // the Bio package (which owns the merged lib) may be absent — sets then load with
  // structurally-complete placeholders whose monomer links resolve lazily elsewhere
  private async tryGetMonomerLib(): Promise<IMonomerLib> {
    try {
      return (await getMonomerLibHelper()).getMonomerLib();
    } catch {
      return {getMonomer: () => null} as unknown as IMonomerLib;
    }
  }
}

function normName(name: string): string {
  return name.endsWith('.json') ? name.substring(0, name.length - 5) : name;
}

function sameList(a: string[], b: string[]): boolean {
  return a.length === b.length && a.every((v, i) => v === b[i]);
}

function monomerKey(polymerType: string, symbol: string): string {
  return `${polymerType} ${symbol}`;
}

function checkBatchReport(report: DG.DomainBatchReport, what: string): void {
  if (report.error != null || report.errorCount > 0)
    throw new Error(`Failed to save ${what}: ${report.error ?? `${report.errorCount} row(s) failed validation`}`);
}

function parseJson<T>(value: string | null | undefined, fallback: T): T {
  if (value == null || value === '')
    return fallback;
  return JSON.parse(value) as T;
}

function rowToMonomer(r: MonomerRow): Monomer {
  const m: Monomer = {
    symbol: r.symbol,
    name: r.name ?? '',
    molfile: r.molfile ?? '',
    author: r.author ?? '',
    id: r.helm_id ?? 0,
    rgroups: parseJson(r.rgroups, []),
    smiles: r.smiles ?? '',
    polymerType: r.polymer_type as PolymerType,
    monomerType: (r.monomer_type ?? 'Undefined') as MonomerType,
    createDate: r.create_date ?? null,
  };
  if (r.natural_analog)
    m.naturalAnalog = r.natural_analog;
  const meta = parseJson<{[k: string]: any} | null>(r.meta, null);
  if (meta != null)
    m.meta = meta;
  return m;
}

function monomersToDf(libraryId: string, monomers: Monomer[]): DG.DataFrame {
  const s = (name: keyof MonomerInsert, get: (m: Monomer) => string | null) =>
    DG.Column.fromList('string', name as string, monomers.map(get));
  return DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'library_id', monomers.map(() => libraryId)),
    s('symbol', (m) => m.symbol),
    s('name', (m) => m.name ?? null),
    s('polymer_type', (m) => m.polymerType),
    s('monomer_type', (m) => m.monomerType ?? null),
    s('molfile', (m) => m.molfile ?? null),
    s('smiles', (m) => m.smiles ?? null),
    s('natural_analog', (m) => m.naturalAnalog ?? null),
    s('author', (m) => m.author ?? null),
    s('create_date', (m) => m.createDate ?? null),
    DG.Column.fromList('int', 'helm_id', monomers.map((m) => typeof m.id === 'number' ? m.id : null)),
    s('rgroups', (m) => JSON.stringify(m.rgroups ?? [])),
    s('meta', (m) => m.meta != null ? JSON.stringify(m.meta) : null),
  ]);
}
