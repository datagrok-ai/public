import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {funcs} from '../package-api';
import {Vault} from '../cdd-vault-api';
import {BATCHES_TAB, CDDVaultSearchType, MOLECULES_TAB, PROTOCOLS_TAB} from '../constants';

// 8374 - id of test Datagrok vault
const TEST_VAULT_ID = 8374;
const TEST_STRUCTURE = 'c1ccccc1';

const APP_NODE_PATH = ['Apps', 'Chem', 'CDD Vault'];

function getCddNode(): DG.TreeViewGroup {
  let node = grok.shell.browsePanel.mainTree;
  for (const name of APP_NODE_PATH)
    node = node.getOrCreateGroup(name);
  return node;
}

category('CDDVaultLink: app', () => {
  test('app launches and shows initial statistics', async () => {
    const view = await funcs.cddVaultApp();
    expect(view != null, true, 'cddVaultApp returned no view');
    // initial stats are rendered as a table inside the app view
    await awaitCheck(() => view.root.getElementsByTagName('table').length > 0,
      'Initial statistics have not been rendered', 30000);
  }, {timeout: 60000});

  test('open Molecules tab', async () => {
    const node = getCddNode();
    node.expanded = true;
    // wait for the vaults to be loaded into the tree
    await awaitCheck(() => node.items.length > 0,
      'CDD Vault tree has not been populated', 30000);
    const vaultNode = node.items[0] as DG.TreeViewGroup;
    vaultNode.expanded = true;
    await awaitCheck(() => vaultNode.items.find((it) => it.text === MOLECULES_TAB) !== undefined,
      `${MOLECULES_TAB} node has not been created`, 10000);
    vaultNode.items.find((it) => it.text === MOLECULES_TAB)!.root.click();
    await awaitCheck(() => (grok.shell.tv?.dataFrame?.rowCount ?? 0) > 0,
      'Molecules dataframe has not been loaded', 30000);
  }, {timeout: 60000});

  test('open Batches tab', async () => {
    const node = getCddNode();
    node.expanded = true;
    await awaitCheck(() => node.items.length > 0,
      'CDD Vault tree has not been populated', 30000);
    const vaultNode = node.items[0] as DG.TreeViewGroup;
    vaultNode.expanded = true;
    await awaitCheck(() => vaultNode.items.find((it) => it.text === BATCHES_TAB) !== undefined,
      `${BATCHES_TAB} node has not been created`, 10000);
    vaultNode.items.find((it) => it.text === BATCHES_TAB)!.root.click();
    await awaitCheck(() => (grok.shell.tv?.dataFrame?.rowCount ?? 0) > 0,
      'Batches dataframe has not been loaded', 30000);
  }, {timeout: 60000});
});

category('CDDVaultLink: functions', () => {
  let vaultId = TEST_VAULT_ID;

  before(async () => {
    // resolve the first available vault id at runtime so the tests don't depend
    // on a hardcoded id from a specific CDD account
    const vaults = JSON.parse(await funcs.getVaults()) as Vault[];
    if (vaults.length > 0)
      vaultId = vaults[0].id;
  });

  test('getVaults returns at least one vault', async () => {
    const vaults = JSON.parse(await funcs.getVaults()) as Vault[];
    expect(Array.isArray(vaults), true, 'getVaults did not return an array');
    expect(vaults.length > 0, true, 'getVaults returned an empty list');
  }, {timeout: 60000});

  test('getMolecules preview returns a dataframe', async () => {
    const df = await funcs.getMolecules(vaultId, '');
    expect(df != null, true, 'getMolecules returned null');
    expect(df.rowCount > 0, true, 'getMolecules returned an empty dataframe');
    const expectedColumns = ['id', 'name', 'smiles', 'class', 'created_at', 'modified_at', 'synonyms',
      'registration_type', 'registration_form', 'projects', 'owner', 'cxsmiles', 'inchi', 'inchi_key',
      'iupac_name', 'molfile', 'molecular_weight', 'log_p', 'log_d', 'log_s', 'num_aromatic_rings',
      'num_h_bond_donors', 'num_h_bond_acceptors', 'num_rule_of_5_violations', 'formula', 'isotope_formula',
      'p_k_a', 'p_k_a_type', 'p_k_a_basic', 'exact_mass', 'heavy_atom_count', 'composition',
      'isotope_composition', 'topological_polar_surface_area', 'num_rotatable_bonds', 'cns_mpo_score',
      'bbb2_score', 'fsp3', 'molecule_batch_identifiers', 'Stage', 'MW'];
    const actualColumns = df.columns.names();
    for (const col of expectedColumns)
      expect(actualColumns.includes(col), true, `getMolecules dataframe is missing column '${col}'`);
  }, {timeout: 120000});

  test('getBatches preview returns a dataframe', async () => {
    const df = await funcs.getBatches(vaultId);
    expect(df != null, true, 'getBatches returned null');
    expect(df.rowCount > 0, true, 'getBatches returned an empty dataframe');
    const expectedColumns = ['id', 'name', 'class', 'created_at', 'modified_at', 'molecule_batch_identifier', 'owner',
      'projects', 'salt_name', 'formula_weight', 'molecule_id', 'molecule_smiles'];
    const actualColumns = df.columns.names();
    for (const col of expectedColumns)
      expect(actualColumns.includes(col), true, `getBatches dataframe is missing column '${col}'`);
  }, {timeout: 120000});

  test('getProtocolsAsync returns a non-empty list', async () => {
    const json = await funcs.getProtocolsAsync(vaultId, 5);
    expect(json !== '' && json != null, true, 'getProtocolsAsync returned an empty string');
    const protocols = JSON.parse(json) as {id: number}[];
    expect(Array.isArray(protocols), true, 'getProtocolsAsync did not return an array');
    expect(protocols.length > 0, true, 'getProtocolsAsync returned an empty list');
    const expectedProtocolIds = [110458, 110475, 110477, 110481, 110484];
    const actualIds = protocols.map((p) => p.id);
    for (const id of expectedProtocolIds)
      expect(actualIds.includes(id), true, `getProtocolsAsync is missing protocol id ${id}`);
  }, {timeout: 5 * 60000 + 30000});

  test('getCollectionsAsync returns a list', async () => {
    const json = await funcs.getCollectionsAsync(vaultId, 5);
    expect(json !== '' && json != null, true, 'getCollectionsAsync returned an empty string');
    const collections = JSON.parse(json) as {id: number}[];
    expect(Array.isArray(collections), true, 'getCollectionsAsync did not return an array');
    const expectedCollectionIds = [780976];
    const actualIds = collections.map((c) => c.id);
    for (const id of expectedCollectionIds)
      expect(actualIds.includes(id), true, `getCollectionsAsync is missing collection id ${id}`);
  }, {timeout: 5 * 60000 + 30000});

  test('getSavedSearches returns a list', async () => {
    const json = await funcs.getSavedSearches(vaultId);
    expect(json !== '' && json != null, true, 'getSavedSearches returned an empty string');
    const searches = JSON.parse(json) as {id: number}[];
    expect(Array.isArray(searches), true, 'getSavedSearches did not return an array');
    const expectedSearchIds = [19013281];
    const actualIds = searches.map((s) => s.id);
    for (const id of expectedSearchIds)
      expect(actualIds.includes(id), true, `getSavedSearches is missing saved search id ${id}`);
  }, {timeout: 60000});

  test('cDDVaultSearch substructure preview returns a dataframe', async () => {
    const df = await funcs.cDDVaultSearchAsync(vaultId, TEST_STRUCTURE,
      CDDVaultSearchType.SUBSTRUCTURE, 0, null, null);
    expect(df != null, true, 'cDDVaultSearch returned null');
    expect(df.rowCount, 992, `cDDVaultSearch returned ${df.rowCount} rows, expected 992`);
  }, {timeout: 120000});
});
