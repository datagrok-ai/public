import * as DG from 'datagrok-api/dg';

import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  IMonomerLib, IMonomerLibProvider, Monomer, MonomerSet, MonomerSetPlaceholder,
} from '@datagrok-libraries/bio/src/types/monomer-library';
import {PolymerType} from '@datagrok-libraries/bio/src/helm/types';

import {DOMAIN_DB_PROVIDER_NAME, DomainDbLibraryProvider} from '../domain-library-provider';
import {monomersDb} from '../generated/db';

const runId = Math.random().toString(36).substring(2, 8);
const libName = (tag: string) => `mddb-test-${runId}-${tag}`;

const TEST_MOLFILE = `
  Mrv1810 01011900002D

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 N   0  0
    0.8250    0.0000    0.0000 C   0  0
    1.6500    0.0000    0.0000 O   0  0
  1  2  1  0
  2  3  1  0
M  END`;

function mkMonomer(symbol: string, polymerType: PolymerType = 'PEPTIDE', extra: Partial<Monomer> = {}): Monomer {
  return {
    symbol,
    name: `Test monomer ${symbol}`,
    molfile: TEST_MOLFILE,
    author: 'monomerDomainDB tests',
    id: 0,
    rgroups: [
      {capGroupSmiles: '[*:1][H]', alternateId: 'R1-H', capGroupName: 'H', label: 'R1'},
      {capGroupSmiles: 'O[*:2]', alternateId: 'R2-OH', capGroupName: 'OH', label: 'R2'},
    ],
    smiles: 'C(C(=O)O)N',
    polymerType: polymerType,
    monomerType: 'Backbone',
    createDate: '2026-09-01T00:00:00.000Z',
    ...extra,
  };
}

function checkMonomer(actual: Monomer | undefined, expected: Monomer): void {
  if (actual == null)
    throw new Error(`Monomer '${expected.symbol}' not found`);
  expect(actual.symbol, expected.symbol);
  expect(actual.name, expected.name);
  expect(actual.molfile, expected.molfile);
  expect(actual.smiles, expected.smiles);
  expect(actual.polymerType, expected.polymerType);
  expect(actual.monomerType, expected.monomerType);
  expect(actual.author, expected.author);
  expect(actual.createDate, expected.createDate);
  expect(actual.id, expected.id);
  expect(JSON.stringify(actual.rgroups), JSON.stringify(expected.rgroups));
  expect(JSON.stringify(actual.meta ?? null), JSON.stringify(expected.meta ?? null));
  expect(actual.naturalAnalog ?? null, expected.naturalAnalog ?? null);
}

async function findLibraryId(name: string): Promise<string | null> {
  const row = await monomersDb.libraries.first({filter: DG.cond('name', '=', name)});
  return row?.id ?? null;
}

category('monomerDomainDB', () => {
  const provider: DomainDbLibraryProvider = DomainDbLibraryProvider.getInstance();

  after(async () => {
    for (const client of [monomersDb.libraries, monomersDb.monomerSets]) {
      for (;;) {
        const report = await client.deleteWhere(DG.cond('name', 'like', 'mddb-test-%'));
        if (!report.hasMore)
          break;
      }
    }
  });

  test('provider is discoverable by role', async () => {
    const funcs = DG.Func.find({meta: {role: 'monomer-lib-provider'}});
    const ours = funcs.find((f) => f.name === 'getMonomerDomainDbProvider');
    expect(ours != null, true);
    const discovered = await ours!.apply({}) as IMonomerLibProvider;
    expect(discovered.name, DOMAIN_DB_PROVIDER_NAME);
  });

  test('library round-trip preserves all monomer fields', async () => {
    const name = libName('roundtrip');
    const full = mkMonomer('Full', 'PEPTIDE', {
      naturalAnalog: 'A',
      meta: {colors: {default: {line: '#123456', text: '#654321', background: '#ffffff'}}, note: 'extra'},
      id: 42,
    });
    const rna = mkMonomer('rA', 'RNA', {monomerType: 'Branch'});
    try {
      await provider.addOrUpdateLibrary(name, [full, rna]);
      expect((await provider.listLibraries()).includes(name), true);
      const data = await provider.getSingleLibrary(name);
      checkMonomer(data?.['PEPTIDE']?.['Full'], full);
      checkMonomer(data?.['RNA']?.['rA'], rna);
    } finally {
      await provider.deleteLibrary(name).catch(() => {});
    }
  });

  test('getLibraryAsString round-trips through addOrUpdateLibraryString', async () => {
    const src = libName('str-src');
    const dst = libName('str-dst');
    const monomers = [mkMonomer('S1'), mkMonomer('S2', 'CHEM', {monomerType: 'Undefined'})];
    try {
      await provider.addOrUpdateLibrary(src, monomers);
      const content = await provider.getLibraryAsString(src);
      const parsed = JSON.parse(content);
      expect(Array.isArray(parsed), true);
      expect(parsed.length, 2);
      await provider.addOrUpdateLibraryString(dst, content);
      const data = await provider.getSingleLibrary(dst);
      checkMonomer(data?.['PEPTIDE']?.['S1'], monomers[0]);
      checkMonomer(data?.['CHEM']?.['S2'], monomers[1]);
    } finally {
      await provider.deleteLibrary(src).catch(() => {});
      await provider.deleteLibrary(dst).catch(() => {});
    }
  });

  test('addOrUpdateLibrary replaces the full monomer set', async () => {
    const name = libName('replace');
    try {
      await provider.addOrUpdateLibrary(name, [mkMonomer('A'), mkMonomer('B')]);
      await provider.addOrUpdateLibrary(name, [mkMonomer('B', 'PEPTIDE', {name: 'Updated B'}), mkMonomer('C')]);
      const data = await provider.getSingleLibrary(name);
      const symbols = Object.keys(data?.['PEPTIDE'] ?? {}).sort();
      expect(JSON.stringify(symbols), JSON.stringify(['B', 'C']));
      expect(data?.['PEPTIDE']?.['B']?.name, 'Updated B');
    } finally {
      await provider.deleteLibrary(name).catch(() => {});
    }
  });

  test('updateOrAddMonomersInLibrary merges', async () => {
    const name = libName('merge');
    try {
      await provider.addOrUpdateLibrary(name, [mkMonomer('A'), mkMonomer('B')]);
      await provider.updateOrAddMonomersInLibrary(name,
        [mkMonomer('B', 'PEPTIDE', {smiles: 'CCO'}), mkMonomer('D')]);
      const data = await provider.getSingleLibrary(name);
      const symbols = Object.keys(data?.['PEPTIDE'] ?? {}).sort();
      expect(JSON.stringify(symbols), JSON.stringify(['A', 'B', 'D']));
      expect(data?.['PEPTIDE']?.['B']?.smiles, 'CCO');
    } finally {
      await provider.deleteLibrary(name).catch(() => {});
    }
  });

  test('deleteMonomersFromLibrary removes only the named monomers', async () => {
    const name = libName('del-monomers');
    try {
      await provider.addOrUpdateLibrary(name,
        [mkMonomer('A'), mkMonomer('B'), mkMonomer('A', 'RNA', {monomerType: 'Branch'})]);
      await provider.deleteMonomersFromLibrary(name, [{polymerType: 'PEPTIDE', symbol: 'A'}]);
      const data = await provider.getSingleLibrary(name);
      expect(data?.['PEPTIDE']?.['A'] == null, true);
      expect(data?.['PEPTIDE']?.['B'] != null, true);
      expect(data?.['RNA']?.['A'] != null, true);
    } finally {
      await provider.deleteLibrary(name).catch(() => {});
    }
  });

  test('repeated addOrUpdateLibrary is idempotent', async () => {
    const name = libName('idempotent');
    const monomers = [mkMonomer('A'), mkMonomer('B')];
    try {
      await provider.addOrUpdateLibrary(name, monomers);
      await provider.addOrUpdateLibrary(name, monomers);
      const libId = await findLibraryId(name);
      expect(libId != null, true);
      expect(await monomersDb.monomers.count(DG.cond('library_id', '=', libId!)), 2);
    } finally {
      await provider.deleteLibrary(name).catch(() => {});
    }
  });

  test('deleteLibrary cascades to its monomers', async () => {
    const name = libName('cascade');
    try {
      await provider.addOrUpdateLibrary(name, [mkMonomer('A'), mkMonomer('B')]);
      const libId = await findLibraryId(name);
      expect(libId != null, true);
      await provider.deleteLibrary(name);
      expect((await provider.listLibraries()).includes(name), false);
      expect(await findLibraryId(name), null);
      expect(await monomersDb.monomers.count(DG.cond('library_id', '=', libId!)), 0);
    } finally {
      await provider.deleteLibrary(name).catch(() => {});
    }
  });

  test('monomer sets round-trip', async () => {
    const name = libName('set');
    const stubLib = {getMonomer: () => null} as unknown as IMonomerLib;
    const links = [{source: 'HELMCoreLibrary.json', symbol: 'A', notes: 'test link'}];
    const set = new MonomerSet('test set description',
      [new MonomerSetPlaceholder(stubLib, 'X', 'PEPTIDE', 'Backbone', links)]);
    try {
      await provider.addOrUpdateSet(name, set);
      expect((await provider.listSets()).includes(name), true);
      const [loaded] = await provider.loadSets([name]);
      expect(loaded.description, 'test set description');
      expect(loaded.placeholders.length, 1);
      expect(loaded.placeholders[0].symbol, 'X');
      expect(loaded.placeholders[0].polymerType, 'PEPTIDE');
      expect(loaded.placeholders[0].monomerType, 'Backbone');
      expect(JSON.stringify(loaded.placeholders[0].monomerLinks), JSON.stringify(links));
    } finally {
      await provider.deleteSet(name).catch(() => {});
    }
  });

  test('monomer set string round-trip and delete', async () => {
    const name = libName('set-str');
    const content = JSON.stringify({
      description: 'from string',
      placeholders: {
        'Y': {polymerType: 'RNA', monomerType: 'Branch', set: [{source: 'lib.json', symbol: 'rA', notes: ''}]},
      },
    });
    try {
      await provider.addOrUpdateSetString(name, content);
      const [loaded] = await provider.loadSets([name]);
      expect(loaded.description, 'from string');
      expect(loaded.placeholders[0].symbol, 'Y');
      expect(loaded.placeholders[0].polymerType, 'RNA');
      await provider.deleteSet(name);
      expect((await provider.listSets()).includes(name), false);
    } finally {
      await provider.deleteSet(name).catch(() => {});
    }
  });

  test('invalid monomer rejects the write and preserves existing content', async () => {
    const name = libName('invalid');
    try {
      await provider.addOrUpdateLibrary(name, [mkMonomer('A'), mkMonomer('B')]);
      let error: Error | null = null;
      await provider.addOrUpdateLibrary(name,
        [mkMonomer('C', 'NOT_A_POLYMER_TYPE' as PolymerType)]).catch((e) => { error = e; });
      expect(error != null, true);
      const data = await provider.getSingleLibrary(name);
      const symbols = Object.keys(data?.['PEPTIDE'] ?? {}).sort();
      expect(JSON.stringify(symbols), JSON.stringify(['A', 'B']));
    } finally {
      await provider.deleteLibrary(name).catch(() => {});
    }
  });

  test('concurrent writes are serialized', async () => {
    const name = libName('concurrent');
    try {
      await Promise.all([
        provider.addOrUpdateLibrary(name, [mkMonomer('A'), mkMonomer('B')]),
        provider.addOrUpdateLibrary(name, [mkMonomer('C')]),
      ]);
      const data = await provider.getSingleLibrary(name);
      const symbols = Object.keys(data?.['PEPTIDE'] ?? {}).sort();
      expect(JSON.stringify(symbols), JSON.stringify(['C']));
    } finally {
      await provider.deleteLibrary(name).catch(() => {});
    }
  });

  test('refreshLists picks up external changes and fires onChanged', async () => {
    const name = libName('external');
    let fired = false;
    const sub = provider.onChanged.subscribe(() => { fired = true; });
    try {
      await provider.listLibraries();
      fired = false;
      await monomersDb.libraries.insert({name});
      await provider.refreshLists();
      expect((await provider.listLibraries()).includes(name), true);
      expect(fired, true);
    } finally {
      sub.unsubscribe();
      await provider.deleteLibrary(name).catch(() => {});
    }
  });
}, {owner: 'drizhinashvili@datagrok.ai'});
