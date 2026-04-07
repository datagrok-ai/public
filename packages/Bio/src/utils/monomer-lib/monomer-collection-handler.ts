/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MonomerCollection, IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HelmType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {MonomerSelectionWidget} from '@datagrok-libraries/bio/src/utils/monomer-selection-dialog';

import {MonomerLibManager} from './lib-manager';
import {SEM_TYPES} from '../constants';
import {_package} from '../../package';

/** Wrapper that pairs a collection name with its loaded data. */
export interface MonomerCollectionInfo {
  name: string; // file name (with .json)
  displayName: string; // name without .json
  data: MonomerCollection;
}

const HANDLER_TYPE = 'MonomerCollection';

/**
 * ObjectHandler for one or more monomer collections.
 * When a single collection is selected it renders details in the context panel.
 * When multiple collections are selected it renders merge / bulk-delete actions.
 */
export class MonomerCollectionHandler extends DG.ObjectHandler<MonomerCollectionInfo | MonomerCollectionInfo[]> {
  get type() { return HANDLER_TYPE; }
  get name() { return 'Monomer Collection'; }

  isApplicable(x: any): boolean {
    if (Array.isArray(x))
      return x.length > 0 && x.every((item: any) => this.isSingleCollection(item));
    return this.isSingleCollection(x);
  }

  private isSingleCollection(x: any): boolean {
    return x != null && typeof x === 'object' && 'name' in x && 'displayName' in x && 'data' in x &&
      typeof (x as MonomerCollectionInfo).data?.monomerSymbols !== 'undefined';
  }

  getCaption(x: MonomerCollectionInfo | MonomerCollectionInfo[]): string {
    if (Array.isArray(x))
      return `${x.length} collections`;
    return x.displayName;
  }

  renderProperties(x: MonomerCollectionInfo | MonomerCollectionInfo[], context: any = null): HTMLElement {
    if (Array.isArray(x))
      return this.renderMultiProperties(x);
    return this.renderSingleProperties(x);
  }

  // -- Single collection panel --

  private renderSingleProperties(info: MonomerCollectionInfo): HTMLElement {
    const acc = ui.accordion('Monomer Collection');
    const currentUser = DG.User.current().login;
    const isOwner = info.data.updatedBy === currentUser;

    // Details pane
    acc.addPane('Details', () => {
      const map: {[key: string]: any} = {
        'Name': info.displayName,
      };
      if (info.data.description)
        map['Description'] = info.data.description;
      if (info.data.tags && info.data.tags.length > 0)
        map['Tags'] = info.data.tags.join(', ');
      map['Monomers'] = `${info.data.monomerSymbols.length}`;
      if (info.data.updatedBy)
        map['Updated by'] = info.data.updatedBy;
      if (info.data.updatedOn) {
        try { map['Updated on'] = new Date(info.data.updatedOn).toLocaleString(); } catch { /* ignore */ }
      }
      return ui.tableFromMap(map);
    }, true);

    // Monomers pane
    acc.addPane('Monomers', () => {
      const container = ui.div([], {style: {display: 'flex', flexWrap: 'wrap', gap: '4px'}});
      const monomerLib = this.getMonomerLib();
      for (const symbol of info.data.monomerSymbols) {
        const tag = ui.div([symbol], {
          style: {
            display: 'inline-block', padding: '2px 8px', borderRadius: '4px',
            background: 'var(--grey-1)', border: '1px solid var(--grey-2)',
            fontSize: '11px', color: 'var(--grey-5)', cursor: 'pointer',
          },
        });
        tag.addEventListener('mouseenter', () => {
          if (monomerLib) {
            const tooltipEl = MonomerCollectionHandler.getMonomerTooltipSafe(monomerLib, symbol);
            if (tooltipEl) {
              const rect = tag.getBoundingClientRect();
              ui.tooltip.show(tooltipEl, rect.left, rect.bottom + 4);
            }
          }
        });
        tag.addEventListener('mouseleave', () => ui.tooltip.hide());
        tag.addEventListener('click', () => {
          grok.shell.o = DG.SemanticValue.fromValueType(symbol, SEM_TYPES.MONOMER);
        });
        container.appendChild(tag);
      }
      return container;
    }, true);

    // Actions pane
    acc.addPane('Actions', () => {
      const actions: HTMLElement[] = [];

      const editBtn = ui.button('Edit', async () => {
        if (isOwner) await this.editCollection(info);
      }, isOwner ? 'Edit this collection' : 'Only the author can edit');
      if (!isOwner) MonomerCollectionHandler.disableButton(editBtn);
      actions.push(editBtn);

      const deleteBtn = ui.button('Delete', async () => {
        if (isOwner) await this.deleteCollection(info);
      }, isOwner ? 'Delete this collection' : 'Only the author can delete');
      if (!isOwner) MonomerCollectionHandler.disableButton(deleteBtn);
      actions.push(deleteBtn);

      actions.push(ui.button('Duplicate', async () => {
        await this.duplicateCollection(info);
      }, 'Create a copy of this collection'));

      return ui.divH(actions);
    }, true);

    return acc.root;
  }

  // -- Multi-collection panel --

  private renderMultiProperties(items: MonomerCollectionInfo[]): HTMLElement {
    const acc = ui.accordion('Monomer Collections');
    const currentUser = DG.User.current().login;
    const allOwned = items.every((c) => c.data.updatedBy === currentUser);

    // Summary pane
    acc.addPane('Selection', () => {
      const list = ui.list(items.map((c) => c.displayName));
      return ui.divV([
        ui.divText(`${items.length} collections selected`),
        list,
      ]);
    }, true);

    // Actions pane
    acc.addPane('Actions', () => {
      const actions: HTMLElement[] = [];

      actions.push(ui.button('Merge', async () => {
        await this.mergeCollections(items);
      }, 'Merge selected collections into one'));

      const deleteAllBtn = ui.button('Delete All', async () => {
        if (allOwned) await this.deleteMultipleCollections(items);
      }, allOwned ? 'Delete all selected collections' : 'Only the author can delete (you must be author of all)');
      if (!allOwned) MonomerCollectionHandler.disableButton(deleteAllBtn);
      actions.push(deleteAllBtn);

      return ui.divH(actions);
    }, true);

    return acc.root;
  }

  // -- Actions --

  private async editCollection(info: MonomerCollectionInfo): Promise<void> {
    // Find existing view and trigger edit
    const view = Array.from(grok.shell.views).find((v) => v.name === 'Monomer Collections');
    if (view && (view as any)._mcView) {
      const mcView = (view as any)._mcView as import('./monomer-collections-view').MonomerCollectionsView;
      mcView.editCollectionPublic(info.name, info.data);
    } else {
      // Fallback: open the edit dialog directly
      const libHelper = await MonomerLibManager.getInstance();
      await libHelper.awaitLoaded();
      const monomerLib = libHelper.getMonomerLib();
      this.showEditDialog(info, libHelper, monomerLib!);
    }
  }

  private showEditDialog(info: MonomerCollectionInfo, libHelper: MonomerLibManager, monomerLib: IMonomerLib): void {
    const polymerTypes: PolymerType[] = ['PEPTIDE', 'RNA', 'CHEM'];
    const nameInput = ui.input.string('Name', {value: info.displayName, nullable: false});
    const descInput = ui.input.string('Description', {value: info.data.description ?? '', nullable: true});
    const tagsInput = ui.input.string('Tags', {value: (info.data.tags ?? []).join(', '), nullable: true, placeholder: 'Comma-separated tags'});
    const polymerTypeInput = ui.input.choice('Polymer Type', {items: [...polymerTypes], value: 'PEPTIDE', nullable: false});
    const monomerWidget = new MonomerSelectionWidget(monomerLib, polymerTypeInput.value as PolymerType, [...(info.data.monomerSymbols ?? [])], libHelper);

    // eslint-disable-next-line rxjs/no-ignored-subscription
    polymerTypeInput.onChanged.subscribe(() => monomerWidget.setPolymerType(polymerTypeInput.value as PolymerType));

    const dlg = ui.dialog({title: `Edit Collection: ${info.displayName}`})
      .add(ui.divV([nameInput, descInput, tagsInput, polymerTypeInput, ui.element('hr'), monomerWidget.root]))
      .onOK(async () => {
        const name = nameInput.value?.trim();
        if (!name) { grok.shell.warning('Collection name is required'); return; }
        const selectedMonomers = monomerWidget.getSelectedMonomers();
        if (selectedMonomers.length === 0) { grok.shell.warning('Please select at least one monomer'); return; }
        const tags = this.parseTags(tagsInput.value);
        try {
          if (name !== info.displayName)
            await libHelper.deleteMonomerCollection(info.name);
          await libHelper.addOrUpdateMonomerCollection(name, selectedMonomers, descInput.value ?? undefined, tags.length > 0 ? tags : undefined);
          grok.shell.info(`Collection "${name}" saved successfully`);
          await this.refreshView();
        } catch (err: any) {
          grok.shell.error('Error saving collection');
        }
      })
      .show({resizable: true});
    dlg.root.style.width = '550px';
    monomerWidget.focus();
  }

  private async deleteCollection(info: MonomerCollectionInfo): Promise<void> {
    ui.dialog({title: 'Delete Collection'})
      .add(ui.divText(`Are you sure you want to delete "${info.displayName}"?`))
      .onOK(async () => {
        try {
          const libHelper = await MonomerLibManager.getInstance();
          await libHelper.deleteMonomerCollection(info.name);
          grok.shell.info(`Collection "${info.displayName}" deleted`);
          await this.refreshView();
        } catch (err: any) {
          grok.shell.error('Error deleting collection');
        }
      })
      .show();
  }

  private async duplicateCollection(info: MonomerCollectionInfo): Promise<void> {
    const nameInput = ui.input.string('New Name', {value: `${info.displayName} (copy)`, nullable: false});
    ui.dialog({title: 'Duplicate Collection'})
      .add(nameInput)
      .onOK(async () => {
        const name = nameInput.value?.trim();
        if (!name) { grok.shell.warning('Name is required'); return; }
        try {
          const libHelper = await MonomerLibManager.getInstance();
          await libHelper.addOrUpdateMonomerCollection(name, [...info.data.monomerSymbols], info.data.description, info.data.tags ? [...info.data.tags] : undefined);
          grok.shell.info(`Collection "${name}" created`);
          await this.refreshView();
        } catch (err: any) {
          grok.shell.error('Error duplicating collection');
        }
      })
      .show();
  }

  private async mergeCollections(items: MonomerCollectionInfo[]): Promise<void> {
    const nameInput = ui.input.string('Merged Name', {value: 'Merged Collection', nullable: false});
    const descInput = ui.input.string('Description', {nullable: true, placeholder: 'Optional description'});
    ui.dialog({title: `Merge ${items.length} Collections`})
      .add(ui.divV([
        ui.divText(`Merging: ${items.map((c) => c.displayName).join(', ')}`),
        nameInput,
        descInput,
      ]))
      .onOK(async () => {
        const name = nameInput.value?.trim();
        if (!name) { grok.shell.warning('Name is required'); return; }
        try {
          const allSymbols = [...new Set(items.flatMap((c) => c.data.monomerSymbols))];
          const allTags = [...new Set(items.flatMap((c) => c.data.tags ?? []))];
          const libHelper = await MonomerLibManager.getInstance();
          await libHelper.addOrUpdateMonomerCollection(name, allSymbols, descInput.value ?? undefined, allTags.length > 0 ? allTags : undefined);
          grok.shell.info(`Merged collection "${name}" created with ${allSymbols.length} monomers`);
          await this.refreshView();
        } catch (err: any) {
          grok.shell.error('Error merging collections');
        }
      })
      .show();
  }

  private async deleteMultipleCollections(items: MonomerCollectionInfo[]): Promise<void> {
    ui.dialog({title: 'Delete Collections'})
      .add(ui.divText(`Are you sure you want to delete ${items.length} collections?\n${items.map((c) => c.displayName).join(', ')}`))
      .onOK(async () => {
        try {
          const libHelper = await MonomerLibManager.getInstance();
          for (const item of items)
            await libHelper.deleteMonomerCollection(item.name);
          grok.shell.info(`${items.length} collections deleted`);
          await this.refreshView();
        } catch (err: any) {
          grok.shell.error('Error deleting collections');
        }
      })
      .show();
  }

  private async refreshView(): Promise<void> {
    const view = Array.from(grok.shell.views).find((v) => v.name === 'Monomer Collections');
    if (view && (view as any)._mcView) {
      const mcView = (view as any)._mcView as import('./monomer-collections-view').MonomerCollectionsView;
      await mcView.refresh();
    }
  }

  private parseTags(input: string | null): string[] {
    if (!input) return [];
    return input.split(',').map((t) => t.trim()).filter((t) => t.length > 0);
  }

  /** Gets the current monomer lib (synchronous, returns null if not loaded). */
  private getMonomerLib(): IMonomerLib | null {
    return _package.monomerLib ?? null;
  }

  /** Gets a monomer tooltip element, searching across all HELM types. */
  private static getMonomerTooltipSafe(monomerLib: IMonomerLib, symbol: string): HTMLElement | null {
    const helmTypes: HelmType[] = [HelmTypes.AA, HelmTypes.NUCLEOTIDE, HelmTypes.CHEM, HelmTypes.BLOB];
    for (const ht of helmTypes) {
      try {
        const wem = monomerLib.getWebEditorMonomer(ht, symbol);
        if (wem) return monomerLib.getTooltip(ht, symbol);
      } catch { /* skip */ }
    }
    return null;
  }

  /** Applies custom disabled styling to a button (greyed out, not clickable). */
  private static disableButton(btn: HTMLElement): void {
    btn.style.opacity = '0.4';
    btn.style.cursor = 'default';
  }

  /** Builds a context menu for a single collection. */
  static buildContextMenu(menu: DG.Menu, info: MonomerCollectionInfo, handler: MonomerCollectionHandler): void {
    const currentUser = DG.User.current().login;
    const isOwner = info.data.updatedBy === currentUser;

    menu.item('Edit', () => { if (isOwner) handler.editCollection(info); }, null,
      {isEnabled: () => isOwner ? null : 'Only the author can edit'});
    menu.item('Delete', () => { if (isOwner) handler.deleteCollection(info); }, null,
      {isEnabled: () => isOwner ? null : 'Only the author can delete'});
    menu.item('Duplicate', () => handler.duplicateCollection(info));
  }

  /** Builds a context menu for multiple collections. */
  static buildMultiContextMenu(menu: DG.Menu, items: MonomerCollectionInfo[], handler: MonomerCollectionHandler): void {
    const currentUser = DG.User.current().login;
    const allOwned = items.every((c) => c.data.updatedBy === currentUser);

    menu.item('Merge', () => handler.mergeCollections(items));
    menu.item('Delete All', () => { if (allOwned) handler.deleteMultipleCollections(items); }, null,
      {isEnabled: () => allOwned ? null : 'Only the author can delete (you must be author of all)'});
  }
}
