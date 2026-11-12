import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {MonomerCollection, IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {MonomerSelectionWidget} from '@datagrok-libraries/bio/src/utils/monomer-selection-dialog';
import {ChoiceInput, Dialog, Form, TagsInput, TextInput, computed, div, signal, span} from '@datagrok-libraries/u2';

import {MonomerLibManager} from './lib-manager';

/** A collection file and its content, as a card or the context panel holds it. */
export interface CollectionRef {
  /** The file name, `.json` included. */
  name: string;
  data: MonomerCollection;
}

/** Who changes collections: the user's login, and whether the user may edit the file share the collections are kept
 * in (Edit on `System:AppData`), which lets them change anyone's collection. */
export interface CollectionAccess {
  user: string;
  shareEditor: boolean;
}

const POLYMER_TYPES: PolymerType[] = ['PEPTIDE', 'RNA', 'CHEM'];

/** The collection dialogs, shared by the Monomer Collections view and the context panel. Each takes
 * `done`, run after the server changed. */
export class CollectionDialogs {
  static displayName(fileName: string): string {
    return fileName.replace(/\.json$/i, '');
  }

  /** The current user's {@link CollectionAccess}, asked of the server: read when a view of the collections opens. */
  static async access(): Promise<CollectionAccess> {
    const user = DG.User.current().login;
    try {
      const share = await grok.dapi.connections.filter('namespace = "System:" and shortName = "AppData"').first();
      return {user, shareEditor: share != null && await grok.dapi.permissions.check(share, 'Edit')};
    }
    catch (e) {
      console.warn('Monomer collections: could not read the AppData share\'s permissions', e);
      return {user, shareEditor: false};
    }
  }

  /** The author, anyone for a collection that names none, and an editor of the share may change it. */
  static canEdit(data: MonomerCollection, access: CollectionAccess): boolean {
    return !data.updatedBy || data.updatedBy === access.user || access.shareEditor;
  }

  /** Why the user may not change `data`, or null. */
  static notEditable(data: MonomerCollection, access: CollectionAccess): string | null {
    return CollectionDialogs.canEdit(data, access) ? null : CollectionDialogs.onlyAuthor(data, 'this collection');
  }

  /** "Only ann or an editor of the AppData share can change Foo". */
  static onlyAuthor(data: MonomerCollection, what: string): string {
    return `Only ${data.updatedBy} or an editor of the AppData share can change ${what}`;
  }

  static count(n: number, what: string): string {
    return `${n} ${what}${n === 1 ? '' : 's'}`;
  }

  /** New Monomer Collection, or Edit Collection when `existing` is given. OK stays disabled, saying
   * why, until the collection has a name no other collection has and a monomer. */
  static async edit(libHelper: MonomerLibManager, monomerLib: IMonomerLib, existing: CollectionRef | null,
    done: () => unknown): Promise<void> {
    const taken = await CollectionDialogs._names(libHelper, existing?.name);
    const shownName = existing === null ? '' : CollectionDialogs.displayName(existing.name);
    const title = existing === null ? 'New Monomer Collection' : `Edit Collection: ${shownName}`;
    const dialog = Dialog.create(title, {name: existing === null ? 'newCollection' : 'editCollection'});
    const name = dialog.build(() => CollectionDialogs._nameInput('Name', shownName, taken));
    const description = dialog.build(() => new TextInput({label: 'Description', name: 'description',
      value: existing?.data.description ?? '', placeholder: 'Optional description'}));
    const tags = dialog.build(() => new TagsInput<string>({label: 'Tags', name: 'tags',
      value: [...existing?.data.tags ?? []], allowNew: true, placeholder: 'Type a tag and press Enter'}));
    const symbols = [...existing?.data.monomerSymbols ?? []];
    // an edited collection searches where its monomers are: Sugars among RNA's
    const type = POLYMER_TYPES.find((t) => {
      const known = new Set(monomerLib.getMonomerSymbolsByType(t));
      return symbols.some((s) => known.has(s));
    }) ?? 'PEPTIDE';
    const monomers = new MonomerSelectionWidget(monomerLib, type, symbols, libHelper);
    const polymerType = dialog.build(() => new ChoiceInput({label: 'Search in', name: 'polymerType',
      items: POLYMER_TYPES, value: type, nullable: false, onInput: (t) => monomers.setPolymerType(t as PolymerType),
      tooltipText: 'The polymer type whose monomers the search below offers'}));
    const picked = signal(monomers.getSelectedMonomers().length);
    const changes = monomers.onChanged.subscribe((symbols) => picked.value = symbols.length);

    dialog.add(() => div([new Form({layout: 'normal'}).addAll([name, description, tags, polymerType]).root,
      div([monomers.root], 'monomer-collection-dialog-monomers')], 'monomer-collection-dialog'))
      .okDisabled(computed(() => picked.value === 0 ? 'Select at least one monomer' : null))
      .onOK(async () => {
        // as the name check reads it: "Foo.json" names Foo's own file, which the rename would delete
        const newName = CollectionDialogs.displayName(name.peek().trim());
        const symbols = monomers.getSelectedMonomers();
        const tagList = tags.peek();
        await libHelper.addOrUpdateMonomerCollection(newName, symbols, description.peek() || undefined,
          tagList.length > 0 ? tagList : undefined);
        // renamed: the old file goes once the new one is written, so a failed write loses nothing
        if (existing !== null && newName !== shownName)
          await libHelper.deleteMonomerCollection(existing.name);
        grok.shell.info(`Collection "${newName}" ${existing === null ? 'created' : 'saved'}`);
        await done();
      }, {text: existing === null ? 'Create' : 'Save'});
    dialog.own(() => changes.unsubscribe());
    dialog.show({width: 560});
    if (existing === null)
      name.focus();
    else
      monomers.focus();
  }

  /** Deletes one collection or several, after the user confirms. */
  static delete(libHelper: MonomerLibManager, refs: CollectionRef[], done: () => unknown): void {
    const names = refs.map((r) => CollectionDialogs.displayName(r.name));
    const text = names.length === 1 ? `Delete the collection "${names[0]}"?` :
      `Delete ${names.length} collections: ${names.join(', ')}?`;
    Dialog.create(names.length === 1 ? 'Delete Collection' : 'Delete Collections', {name: 'deleteCollection'})
      .add(span(text))
      .onOK(async () => {
        const present = await libHelper.listMonomerCollections();
        const gone = refs.filter((r) => !present.includes(r.name));
        for (const ref of refs.filter((r) => present.includes(r.name)))
          await libHelper.deleteMonomerCollection(ref.name);
        const deleted = refs.length - gone.length;
        grok.shell.info(gone.length === refs.length ? 'Already deleted' :
          deleted === 1 && refs.length === 1 ? `Collection "${names[0]}" deleted` :
            `${CollectionDialogs.count(deleted, 'collection')} deleted` +
            (gone.length > 0 ? `, ${gone.length} already gone` : ''));
        await done();
      }, {text: 'Delete'})
      .show();
  }

  /** A copy of a collection under a new name. */
  static async duplicate(libHelper: MonomerLibManager, ref: CollectionRef, done: () => unknown): Promise<void> {
    const taken = await CollectionDialogs._names(libHelper);
    const dialog = Dialog.create('Duplicate Collection', {name: 'duplicateCollection'});
    const name = dialog.build(() => CollectionDialogs._nameInput('New Name',
      `${CollectionDialogs.displayName(ref.name)} (copy)`, taken));
    dialog.add(() => new Form({layout: 'normal'}).add(name))
      .onOK(async () => {
        const newName = CollectionDialogs.displayName(name.peek().trim());
        await libHelper.addOrUpdateMonomerCollection(newName, [...ref.data.monomerSymbols], ref.data.description,
          ref.data.tags ? [...ref.data.tags] : undefined);
        grok.shell.info(`Collection "${newName}" created`);
        await done();
      }, {text: 'Duplicate'})
      .show();
  }

  /** One collection with every monomer and tag of `refs`. */
  static async merge(libHelper: MonomerLibManager, refs: CollectionRef[], done: () => unknown): Promise<void> {
    const taken = await CollectionDialogs._names(libHelper);
    const dialog = Dialog.create(`Merge ${refs.length} Collections`, {name: 'mergeCollections'});
    const name = dialog.build(() => CollectionDialogs._nameInput('Merged Name', 'Merged Collection', taken));
    const description = dialog.build(() =>
      new TextInput({label: 'Description', name: 'description', placeholder: 'Optional description'}));
    const symbols = [...new Set(refs.flatMap((r) => r.data.monomerSymbols))];
    const tags = [...new Set(refs.flatMap((r) => r.data.tags ?? []))];
    dialog.add(() => div([
        span(`${refs.map((r) => CollectionDialogs.displayName(r.name)).join(', ')} — ` +
          CollectionDialogs.count(symbols.length, 'monomer'), 'monomer-collection-dialog-note'),
        new Form({layout: 'normal'}).addAll([name, description]).root,
      ], 'monomer-collection-dialog'))
      .onOK(async () => {
        const newName = CollectionDialogs.displayName(name.peek().trim());
        await libHelper.addOrUpdateMonomerCollection(newName, symbols, description.peek() || undefined,
          tags.length > 0 ? tags : undefined);
        grok.shell.info(`Merged collection "${newName}" created with ` +
          CollectionDialogs.count(symbols.length, 'monomer'));
        await done();
      }, {text: 'Merge'})
      .show();
  }

  /** The collections' names as typed, lower-cased, but `except`'s — an edit keeping its own name. */
  private static async _names(libHelper: MonomerLibManager, except?: string): Promise<Set<string>> {
    return new Set((await libHelper.listMonomerCollections()).filter((n) => n !== except)
      .map((n) => CollectionDialogs.displayName(n).toLowerCase()));
  }

  /** A required name that no other collection has: saving would overwrite that one. */
  private static _nameInput(label: string, value: string, taken: Set<string>): TextInput {
    const input = new TextInput({label, name: 'name', required: true, value, placeholder: 'Enter collection name'});
    input.addValidator((v) => taken.has(CollectionDialogs.displayName(v.trim()).toLowerCase()) ?
      'A collection with this name already exists' : null);
    return input;
  }
}
