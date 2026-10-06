/* eslint-disable max-len */
/**
 * AI view functions for the monomer-management views (Manage Monomer Libraries,
 * Manage Monomers, Monomer Collections).
 *
 * The vocabulary is registered ONCE per session through `grok.functions.register`
 * (real platform Funcs, no package.ts entries). Each view gets an instance
 * `getFunctions` override on the wrapper the shell reports (`DG.toJs(view.dart)` —
 * a `DG.View` built with `new View(...)` is NOT the wrapper `grok.shell.v` later
 * returns, so the override must live on the Dart-cached one), plus an
 * `aiDescription` briefing.
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerCollection} from '@datagrok-libraries/bio/src/types/monomer-library';
import {PolymerType} from '@datagrok-libraries/bio/src/helm/types';

import {MonomerLibManager} from './lib-manager';
import {MonomerManager} from './monomer-manager/monomer-manager';
import {DuplicateMonomerManager} from './monomer-manager/duplicate-monomer-manager';
import {deleteMonomerLibraryByName, mergeActiveLibrariesInto, refreshLibraryControlsForm}
  from './library-file-manager/ui';
import type {MonomerCollectionsView} from './monomer-collections-view';

const AI_NAMESPACE = 'BioMonomerAI';

const _funcs = new Map<string, DG.Func>();

/** Registers one vocabulary Func on first use; later calls return the cached one. */
function aiFunc(name: string, signature: string, description: string, run: (...args: any[]) => any): DG.Func {
  let f = _funcs.get(name);
  if (!f) {
    f = grok.functions.register({signature, run, isAsync: true, namespace: AI_NAMESPACE,
      tags: 'bio-monomer-ai-function', options: {description}});
    _funcs.set(name, f);
  }
  return f;
}

/** The wrapper instance `grok.shell.v` reports for [view] — overrides must live there. */
function shellWrapperOf(view: DG.ViewBase): DG.ViewBase {
  return (view.dart ? DG.toJs(view.dart) : null) ?? view;
}

async function libManager(): Promise<MonomerLibManager> {
  const manager = await MonomerLibManager.getInstance();
  await manager.awaitLoaded();
  return manager;
}

// ─────────────────────── Manage Monomer Libraries ───────────────────────

function libraryManagerFunctions(): DG.Func[] {
  return [
    aiFunc('getMonomerLibrariesState', 'dynamic getMonomerLibrariesState(view view)',
      'The monomer libraries - name, storage provider, whether each is active (its monomers loaded), and the duplicate monomer symbols found across libraries with the preferred source of each. Call this first',
      async () => {
        const manager = await libManager();
        const settings = await getUserLibSettings();
        const perProvider = await manager.getAvailableLibrariesPerProvider();
        const libraries: {name: string, provider: string, active: boolean}[] = [];
        for (const [provider, names] of perProvider) {
          for (const name of names)
            libraries.push({name, provider, active: !settings.exclude.includes(name)});
        }
        const duplicates: {polymerType: string, symbol: string, libraries: string[], preferred: string | null}[] = [];
        for (const [polymerType, bySymbol] of Object.entries(manager.duplicateMonomers)) {
          for (const [symbol, monomers] of Object.entries(bySymbol)) {
            duplicates.push({
              polymerType, symbol,
              libraries: monomers.map((m) => m.lib?.source ?? '').filter(Boolean),
              preferred: settings.duplicateMonomerPreferences?.[polymerType]?.[symbol] ?? null,
            });
          }
        }
        return {libraries, duplicates, allDuplicatesResolved: manager.duplicatesHandled};
      }),

    aiFunc('setMonomerLibraryActive', 'dynamic setMonomerLibraryActive(view view, string library, bool active)',
      'Activate or deactivate a monomer library - active libraries contribute their monomers to the merged set used across the platform',
      async (_view: any, library: string, active: boolean) => {
        const manager = await libManager();
        const names = await manager.getAvaliableLibraryNames();
        if (!names.includes(library))
          return {success: false, error: `No library '${library}' — call getMonomerLibrariesState for names`};
        const settings = await getUserLibSettings();
        if (active)
          settings.exclude = settings.exclude.filter((n) => n !== library);
        else if (!settings.exclude.includes(library))
          settings.exclude.push(library);
        await setUserLibSettings(settings);
        await manager.loadMonomerLib(true);
        manager.notifyLibrarySelectionChanged();
        await refreshLibraryControlsForm();
        return {success: true, library, active: !!active};
      }),

    aiFunc('resolveDuplicateMonomer', 'dynamic resolveDuplicateMonomer(view view, string symbol, string library)',
      'For a monomer symbol duplicated across libraries, choose which library its definition should come from',
      async (_view: any, symbol: string, library: string) => {
        const manager = await libManager();
        const entries = Object.entries(manager.duplicateMonomers)
          .filter(([, bySymbol]) => symbol in bySymbol);
        if (!entries.length)
          return {success: false, error: `'${symbol}' is not a duplicated monomer symbol — call getMonomerLibrariesState`};
        const settings = await getUserLibSettings();
        const applied: string[] = [];
        for (const [polymerType, bySymbol] of entries) {
          const sources = bySymbol[symbol].map((m) => m.lib?.source ?? '');
          if (!sources.includes(library))
            return {success: false, error: `'${symbol}' has no definition in '${library}' — its sources are: ${sources.filter(Boolean).join(', ')}`};
          settings.duplicateMonomerPreferences ??= {};
          settings.duplicateMonomerPreferences[polymerType] ??= {};
          settings.duplicateMonomerPreferences[polymerType][symbol] = library;
          applied.push(polymerType);
        }
        await setUserLibSettings(settings);
        manager.assignDuplicatePreferances(settings);
        await (await DuplicateMonomerManager.getInstance()).refresh();
        return {success: true, symbol, library, polymerTypes: applied};
      }),

    aiFunc('deleteMonomerLibrary', 'dynamic deleteMonomerLibrary(view view, string library)',
      'Permanently delete a monomer library file. ALWAYS confirm with the user first',
      async (_view: any, library: string) => {
        const manager = await libManager();
        const names = await manager.getAvaliableLibraryNames();
        if (!names.includes(library))
          return {success: false, error: `No library '${library}' — call getMonomerLibrariesState for names`};
        await deleteMonomerLibraryByName(library);
        await refreshLibraryControlsForm();
        return {success: true, deleted: library};
      }),

    aiFunc('mergeActiveMonomerLibraries', 'dynamic mergeActiveMonomerLibraries(view view, string newLibraryName, string storage)',
      'Merge the currently active libraries into one new library and switch to it (the previous ones are deactivated). Requires all duplicate symbols to be resolved first. storage is optional - the provider to save into',
      async (_view: any, newLibraryName: string, storage?: string) => {
        const manager = await libManager();
        if (!manager.duplicatesHandled)
          return {success: false, error: 'Duplicate monomer symbols are unresolved — resolve them with resolveDuplicateMonomer first'};
        let fileName = (newLibraryName ?? '').trim();
        if (!fileName)
          return {success: false, error: 'newLibraryName is required'};
        if (!fileName.toLowerCase().endsWith('.json'))
          fileName += '.json';
        const names = await manager.getAvaliableLibraryNames();
        if (names.includes(fileName))
          return {success: false, error: `Library '${fileName}' already exists — pick another name`};
        await mergeActiveLibrariesInto(fileName, storage ?? undefined);
        return {success: true, library: fileName};
      }),

    aiFunc('openMonomerManager', 'dynamic openMonomerManager(view view, string library)',
      'Open the Manage Monomers view to browse and edit the monomers of a library (optional - defaults to the first one)',
      async (_view: any, library?: string) => {
        const manager = await libManager();
        if (library && !(await manager.getAvaliableLibraryNames()).includes(library))
          return {success: false, error: `No library '${library}' — call getMonomerLibrariesState for names`};
        grok.shell.v = await (await MonomerManager.getInstance()).getViewRoot(library ?? undefined);
        return {success: true, note: 'The Manage Monomers view is open'};
      }),
  ];
}

export function attachMonomerLibrariesAi(view: DG.ViewBase): void {
  const w = shellWrapperOf(view);
  // temporary before api release
  (w as any).aiDescription = 'Manage Monomer Libraries — HELM monomer library management: activate/deactivate libraries, ' +
    'upload, delete and merge them, and resolve monomer symbols duplicated across libraries. ' +
    'Act through the view functions (search list_view_functions with "monomer"): getMonomerLibrariesState ' +
    '(call first), setMonomerLibraryActive, resolveDuplicateMonomer, deleteMonomerLibrary (confirm with the ' +
    'user first), mergeActiveMonomerLibraries, openMonomerManager.';
  w.getFunctions = () => libraryManagerFunctions();
}

// ─────────────────────── Manage Monomers ───────────────────────

const managerByView = new WeakMap<object, MonomerManager>();

function monomerManagerOf(view: any): MonomerManager {
  const v = view?.jsView ?? view;
  const mgr = (v != null ? managerByView.get(v) : null) ?? (view != null ? managerByView.get(view) : null);
  if (!mgr)
    throw new Error('The current view is not the Manage Monomers view');
  return mgr;
}

function monomerManagerFunctions(): DG.Func[] {
  return [
    aiFunc('getMonomerManagerState', 'dynamic getMonomerManagerState(view view)',
      'The Manage Monomers view state - the open library, its monomer symbols, and all available libraries. Call this first',
      async (view: any) => {
        const mgr = monomerManagerOf(view);
        const symbols = mgr.monomerSymbols;
        return {
          library: mgr.currentLibrary,
          monomers: symbols.length,
          symbols: symbols.slice(0, 300),
          ...(symbols.length > 300 ? {note: `Showing 300 of ${symbols.length} symbols`} : {}),
          availableLibraries: await (await libManager()).getAvaliableLibraryNames(),
        };
      }),

    aiFunc('openMonomerLibrary', 'dynamic openMonomerLibrary(view view, string library)',
      'Switch the Manage Monomers view to another library',
      async (view: any, library: string) => {
        const mgr = monomerManagerOf(view);
        const names = await (await libManager()).getAvaliableLibraryNames();
        if (!names.includes(library))
          return {success: false, error: `No library '${library}' — the libraries are: ${names.join(', ')}`};
        grok.shell.v = await mgr.getViewRoot(library);
        return {success: true, library};
      }),

    aiFunc('getMonomerDetails', 'dynamic getMonomerDetails(view view, string symbol)',
      'Full detail of one monomer of the open library - structure, types, natural analog and R-groups',
      async (view: any, symbol: string) => {
        const mgr = monomerManagerOf(view);
        const monomer = await mgr.getMonomerBySymbol(symbol);
        if (!monomer)
          return {success: false, error: `No monomer '${symbol}' in the open library — call getMonomerManagerState for symbols`};
        return {
          symbol: monomer.symbol,
          name: monomer.name,
          smiles: monomer.smiles,
          polymerType: monomer.polymerType,
          monomerType: monomer.monomerType,
          ...(monomer.naturalAnalog ? {naturalAnalog: monomer.naturalAnalog} : {}),
          rGroups: (monomer.rgroups ?? []).map((r) => ({label: r.label, capGroupSmiles: r.capGroupSmiles})),
          ...(monomer.author ? {author: monomer.author} : {}),
          ...(monomer.createDate ? {created: monomer.createDate} : {}),
        };
      }),

    aiFunc('editMonomerInForm', 'dynamic editMonomerInForm(view view, string symbol)',
      'Load a monomer of the open library into the edit form so the user (or setMonomerFormFields) can change it',
      async (view: any, symbol: string) => {
        const mgr = monomerManagerOf(view);
        const row = mgr.findMonomerRow(symbol);
        if (row < 0)
          return {success: false, error: `No monomer '${symbol}' in the open library — call getMonomerManagerState for symbols`};
        await mgr.editMonomer(mgr.tableView!.dataFrame.rows.get(row));
        return {success: true, symbol, note: 'The monomer is loaded into the edit form'};
      }),

    aiFunc('setMonomerFormFields', 'dynamic setMonomerFormFields(view view, map fields)',
      'Set fields of the monomer edit form. Keys: molecule (SMILES), symbol, name, polymerType (PEPTIDE, RNA, CHEM), monomerType (Backbone, Branch, Terminal), naturalAnalog, id. Returns the validation state',
      async (view: any, fields: {[key: string]: unknown}) => {
        const mgr = monomerManagerOf(view);
        if (!fields || !Object.keys(fields).length)
          return {success: false, error: 'Pass fields as a map, e.g. molecule, symbol, name'};
        const unknown = mgr.setMonomerFormFields(fields);
        const validationError = mgr.validateMonomerForm();
        return {success: true, ...(unknown.length ? {unknownFields: unknown} : {}),
          valid: !validationError, ...(validationError ? {validationError} : {})};
      }),

    aiFunc('saveMonomerForm', 'dynamic saveMonomerForm(view view)',
      'Save the monomer currently in the edit form into the open library (adds a new monomer or updates the existing one with that symbol). A dialog asks the user to confirm symbol or structure clashes',
      async (view: any) => {
        const mgr = monomerManagerOf(view);
        const validationError = mgr.validateMonomerForm();
        if (validationError)
          return {success: false, error: validationError};
        await mgr.saveMonomerForm();
        return {success: true, note: 'The monomer was saved (or a clash-confirmation dialog was shown to the user)'};
      }),

    aiFunc('deleteMonomers', 'dynamic deleteMonomers(view view, list symbols)',
      'Delete monomers from the open library by symbol. ALWAYS confirm with the user first',
      async (view: any, symbols: string[]) => {
        const mgr = monomerManagerOf(view);
        if (!Array.isArray(symbols) || !symbols.length)
          return {success: false, error: 'Pass symbols as a non-empty list'};
        const res = await mgr.deleteMonomersBySymbols(symbols);
        if (!res.deleted.length)
          return {success: false, error: `No matching monomers: ${res.missing.join(', ')} — call getMonomerManagerState for symbols`};
        return {success: true, deleted: res.deleted, ...(res.missing.length ? {notFound: res.missing} : {})};
      }),

    aiFunc('fixAllMonomers', 'dynamic fixAllMonomers(view view)',
      'Open the confirmation dialog that re-standardizes every monomer of the open library (smiles, molblocks, R-groups, natural analogs) and saves it. The user confirms in the dialog',
      async (view: any) => {
        await monomerManagerOf(view).fixAllMonomers();
        return {success: true, note: 'A confirmation dialog is shown to the user'};
      }),

    aiFunc('createMonomerLibrary', 'dynamic createMonomerLibrary(view view, string libraryName, string storage)',
      'Create a new empty monomer library and open it. storage is optional - the provider to save into',
      async (view: any, libraryName: string, storage?: string) => {
        const mgr = monomerManagerOf(view);
        let name = (libraryName ?? '').trim();
        if (!name)
          return {success: false, error: 'libraryName is required'};
        if (!name.toLowerCase().endsWith('.json'))
          name += '.json';
        const manager = await libManager();
        if ((await manager.getAvaliableLibraryNames()).includes(name))
          return {success: false, error: `Library '${name}' already exists`};
        const providers = await manager.getProviders();
        const provider = storage ? providers.find((p) => p.name === storage) : providers[0];
        if (!provider)
          return {success: false, error: storage ? `No storage provider '${storage}'` : 'No storage providers available'};
        await mgr.createNewMonomerLib(provider.name, name, []);
        return {success: true, library: name};
      }),
  ];
}

export function attachMonomerManagerAi(view: DG.TableView, manager: MonomerManager): void {
  managerByView.set(view, manager);
  // temporary before api release
  (view as any).aiDescription = 'Manage Monomers — the monomer editor of one HELM library: a grid of its monomers plus an ' +
    'edit form (structure sketcher, symbol, name, types, R-groups). Act through the view functions (search ' +
    'list_view_functions with "monomer"): getMonomerManagerState (call first), openMonomerLibrary, ' +
    'getMonomerDetails, editMonomerInForm / setMonomerFormFields / saveMonomerForm to create or change monomers, ' +
    'deleteMonomers (confirm with the user first), fixAllMonomers, createMonomerLibrary. The standard table ' +
    'view commands are available too.';
  view.getFunctions = () => [...monomerManagerFunctions(), ...DG.View.prototype.getFunctions.call(view)];
}

// ─────────────────────── Monomer Collections ───────────────────────

const collectionsViewByWrapper = new WeakMap<object, MonomerCollectionsView>();

function collectionsViewOf(view: any): MonomerCollectionsView {
  const v = view?.jsView ?? view;
  const mcView = (v != null ? collectionsViewByWrapper.get(v) : null) ??
    (view != null ? collectionsViewByWrapper.get(view) : null);
  if (!mcView)
    throw new Error('The current view is not the Monomer Collections view');
  return mcView;
}

function canEditCollection(c: MonomerCollection): boolean {
  return !c.updatedBy || c.updatedBy === DG.User.current().login;
}

/** Which of [symbols] the merged monomer library does not know (any polymer type). */
async function unknownMonomerSymbols(symbols: string[]): Promise<string[]> {
  const lib = (await libManager()).getMonomerLib();
  const polymerTypes = lib.getPolymerTypes();
  return symbols.filter((s) => !polymerTypes.some((pt) => lib.getMonomer(pt as PolymerType, s) != null));
}

function collectionsFunctions(): DG.Func[] {
  return [
    aiFunc('getMonomerCollectionsState', 'dynamic getMonomerCollectionsState(view view)',
      'The saved monomer collections - name, monomer count, tags, author and whether the current user can edit each. Call this first',
      async () => {
        const manager = await libManager();
        const names = await manager.listMonomerCollections();
        const collections = await Promise.all(names.map(async (name) => {
          const c = await manager.readMonomerCollection(name);
          return {
            name: name.replace(/\.json$/i, ''),
            monomers: c.monomerSymbols?.length ?? 0,
            ...(c.tags?.length ? {tags: c.tags} : {}),
            ...(c.updatedBy ? {author: c.updatedBy} : {}),
            canEdit: canEditCollection(c),
          };
        }));
        return {total: collections.length, collections};
      }),

    aiFunc('getMonomerCollectionDetails', 'dynamic getMonomerCollectionDetails(view view, string name)',
      'Full detail of one monomer collection - its monomer symbols, description and tags',
      async (_view: any, name: string) => {
        const manager = await libManager();
        const fileName = await resolveCollectionName(manager, name);
        if (!fileName)
          return {success: false, error: `No collection '${name}' — call getMonomerCollectionsState for names`};
        const c = await manager.readMonomerCollection(fileName);
        return {
          name: fileName.replace(/\.json$/i, ''),
          monomerSymbols: c.monomerSymbols ?? [],
          ...(c.description ? {description: c.description} : {}),
          ...(c.tags?.length ? {tags: c.tags} : {}),
          ...(c.updatedBy ? {author: c.updatedBy} : {}),
          canEdit: canEditCollection(c),
        };
      }),

    aiFunc('createMonomerCollection', 'dynamic createMonomerCollection(view view, string name, list monomerSymbols, string description, list tags)',
      'Create a monomer collection - a named, tagged set of monomer symbols. description and tags are optional',
      async (view: any, name: string, monomerSymbols: string[], description?: string, tags?: string[]) => {
        const manager = await libManager();
        const clean = (name ?? '').trim();
        if (!clean)
          return {success: false, error: 'name is required'};
        if (!Array.isArray(monomerSymbols) || !monomerSymbols.length)
          return {success: false, error: 'Pass monomerSymbols as a non-empty list of symbols'};
        if (await resolveCollectionName(manager, clean))
          return {success: false, error: `Collection '${clean}' already exists — use updateMonomerCollection to change it`};
        const unknown = await unknownMonomerSymbols(monomerSymbols);
        await manager.addOrUpdateMonomerCollection(clean, monomerSymbols,
          description?.trim() || undefined, tags?.length ? tags : undefined);
        await collectionsViewOf(view).refresh();
        return {success: true, name: clean,
          ...(unknown.length ? {note: `Created, but these symbols are not in the active monomer libraries: ${unknown.join(', ')}`} : {})};
      }),

    aiFunc('updateMonomerCollection', 'dynamic updateMonomerCollection(view view, string name, list monomerSymbols, string description, list tags)',
      'Update a monomer collection - each of monomerSymbols, description and tags is optional and replaces the current value when passed',
      async (view: any, name: string, monomerSymbols?: string[], description?: string, tags?: string[]) => {
        const manager = await libManager();
        const fileName = await resolveCollectionName(manager, name);
        if (!fileName)
          return {success: false, error: `No collection '${name}' — call getMonomerCollectionsState for names`};
        const existing = await manager.readMonomerCollection(fileName);
        if (!canEditCollection(existing))
          return {success: false, error: `Collection '${name}' belongs to ${existing.updatedBy} — only its author can edit it`};
        const symbols = Array.isArray(monomerSymbols) && monomerSymbols.length ? monomerSymbols : existing.monomerSymbols;
        const unknown = await unknownMonomerSymbols(symbols);
        await manager.addOrUpdateMonomerCollection(fileName.replace(/\.json$/i, ''), symbols,
          description?.trim() ? description.trim() : existing.description,
          Array.isArray(tags) && tags.length ? tags : existing.tags);
        await collectionsViewOf(view).refresh();
        return {success: true, name: fileName.replace(/\.json$/i, ''),
          ...(unknown.length ? {note: `Saved, but these symbols are not in the active monomer libraries: ${unknown.join(', ')}`} : {})};
      }),

    aiFunc('deleteMonomerCollection', 'dynamic deleteMonomerCollection(view view, string name)',
      'Permanently delete a monomer collection. ALWAYS confirm with the user first',
      async (view: any, name: string) => {
        const manager = await libManager();
        const fileName = await resolveCollectionName(manager, name);
        if (!fileName)
          return {success: false, error: `No collection '${name}' — call getMonomerCollectionsState for names`};
        const existing = await manager.readMonomerCollection(fileName);
        if (!canEditCollection(existing))
          return {success: false, error: `Collection '${name}' belongs to ${existing.updatedBy} — only its author can delete it`};
        await manager.deleteMonomerCollection(fileName);
        await collectionsViewOf(view).refresh();
        return {success: true, deleted: fileName.replace(/\.json$/i, '')};
      }),
  ];
}

/** Stored collection file name for [name] (with or without the .json suffix), or null. */
async function resolveCollectionName(manager: MonomerLibManager, name: string): Promise<string | null> {
  const names = await manager.listMonomerCollections();
  const q = (name ?? '').trim().toLowerCase();
  return names.find((n) => n.toLowerCase() === q || n.replace(/\.json$/i, '').toLowerCase() === q) ?? null;
}

export function attachMonomerCollectionsAi(view: DG.ViewBase, mcView: MonomerCollectionsView): void {
  const w = shellWrapperOf(view);
  collectionsViewByWrapper.set(w, mcView);
  if (w !== view)
    collectionsViewByWrapper.set(view, mcView);
  // temporary before api release
  (w as any).aiDescription = 'Monomer Collections — named, tagged sets of monomer symbols curated from the monomer ' +
    'libraries. Act through the view functions (search list_view_functions with "collection"): ' +
    'getMonomerCollectionsState (call first), getMonomerCollectionDetails, createMonomerCollection, ' +
    'updateMonomerCollection, deleteMonomerCollection (confirm deletions with the user first).';
  w.getFunctions = () => collectionsFunctions();
}
