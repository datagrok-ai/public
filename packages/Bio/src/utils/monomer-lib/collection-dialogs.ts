import * as grok from 'datagrok-api/grok';

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

const POLYMER_TYPES: PolymerType[] = ['PEPTIDE', 'RNA', 'CHEM'];

/** The collection dialogs, shared by the Monomer Collections view and the context panel. Each takes
 * `done`, run after the server changed. */
export class CollectionDialogs {
  static displayName(fileName: string): string {
    return fileName.replace(/\.json$/i, '');
  }

  /** The author, and anyone for a collection that names none, may change it. */
  static canEdit(data: MonomerCollection, user: string): boolean {
    return !data.updatedBy || data.updatedBy === user;
  }

  /** Why `user` may not change `data`, or null. */
  static notEditable(data: MonomerCollection, user: string): string | null {
    return CollectionDialogs.canEdit(data, user) ? null : `Only ${data.updatedBy} can change this collection`;
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
    const name = CollectionDialogs._nameInput('Name', shownName, taken);
    const description = new TextInput({label: 'Description', name: 'description',
      value: existing?.data.description ?? '', placeholder: 'Optional description'});
    const tags = new TagsInput<string>({label: 'Tags', name: 'tags', value: [...existing?.data.tags ?? []],
      allowNew: true, placeholder: 'Type a tag and press Enter'});
    const symbols = [...existing?.data.monomerSymbols ?? []];
    // an edited collection searches where its monomers are: Sugars among RNA's
    const type = POLYMER_TYPES.find((t) => {
      const known = new Set(monomerLib.getMonomerSymbolsByType(t));
      return symbols.some((s) => known.has(s));
    }) ?? 'PEPTIDE';
    const monomers = new MonomerSelectionWidget(monomerLib, type, symbols, libHelper);
    const polymerType = new ChoiceInput({label: 'Search in', name: 'polymerType', items: POLYMER_TYPES,
      value: type, nullable: false, onInput: (t) => monomers.setPolymerType(t as PolymerType),
      tooltipText: 'The polymer type whose monomers the search below offers'});
    const picked = signal(monomers.getSelectedMonomers().length);
    const changes = monomers.onChanged.subscribe((symbols) => picked.value = symbols.length);

    const title = existing === null ? 'New Monomer Collection' : `Edit Collection: ${shownName}`;
    const dialog = Dialog.create(title, {name: existing === null ? 'newCollection' : 'editCollection'});
    dialog.add(div([new Form({layout: 'normal'}).addAll([name, description, tags, polymerType]).root,
      div([monomers.root], 'monomer-collection-dialog-monomers')], 'monomer-collection-dialog'))
      .okDisabled(computed(() => picked.value === 0 ? 'Select at least one monomer' : null))
      .onOK(async () => {
        const newName = name.value.peek().trim();
        const symbols = monomers.getSelectedMonomers();
        const tagList = tags.value.peek();
        await libHelper.addOrUpdateMonomerCollection(newName, symbols, description.value.peek() || undefined,
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
    const name = CollectionDialogs._nameInput('New Name', `${CollectionDialogs.displayName(ref.name)} (copy)`,
      await CollectionDialogs._names(libHelper));
    Dialog.create('Duplicate Collection', {name: 'duplicateCollection'})
      .add(new Form({layout: 'normal'}).add(name))
      .onOK(async () => {
        const newName = name.value.peek().trim();
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
    const name = CollectionDialogs._nameInput('Merged Name', 'Merged Collection', taken);
    const description = new TextInput({label: 'Description', name: 'description', placeholder: 'Optional description'});
    const symbols = [...new Set(refs.flatMap((r) => r.data.monomerSymbols))];
    const tags = [...new Set(refs.flatMap((r) => r.data.tags ?? []))];
    Dialog.create(`Merge ${refs.length} Collections`, {name: 'mergeCollections'})
      .add(div([
        span(`${refs.map((r) => CollectionDialogs.displayName(r.name)).join(', ')} — ` +
          CollectionDialogs.count(symbols.length, 'monomer'), 'monomer-collection-dialog-note'),
        new Form({layout: 'normal'}).addAll([name, description]).root,
      ], 'monomer-collection-dialog'))
      .onOK(async () => {
        const newName = name.value.peek().trim();
        await libHelper.addOrUpdateMonomerCollection(newName, symbols, description.value.peek() || undefined,
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
