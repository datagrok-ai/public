/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MonomerCollection, IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';

import {MonomerLibManager} from './lib-manager';
import {CollectionDialogs} from './collection-dialogs';
import {Control} from '@datagrok-libraries/u2';
import type {MonomerCollectionsView} from './monomer-collections-view';
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
    const isOwner = CollectionDialogs.canEdit(info.data, currentUser);

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
      }, CollectionDialogs.notEditable(info.data, currentUser) ?? 'Edit this collection');
      if (!isOwner) MonomerCollectionHandler.disableButton(editBtn);
      actions.push(editBtn);

      const deleteBtn = ui.button('Delete', async () => {
        if (isOwner) await this.deleteCollection(info);
      }, CollectionDialogs.notEditable(info.data, currentUser) ?? 'Delete this collection');
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
    const allOwned = items.every((c) => CollectionDialogs.canEdit(c.data, currentUser));

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

      const deleteAllBtn = ui.button(`Delete ${items.length} Collections`, async () => {
        if (allOwned) await this.deleteMultipleCollections(items);
      }, MonomerCollectionHandler._notAllEditable(items, currentUser) ?? 'Delete the selected collections');
      if (!allOwned) MonomerCollectionHandler.disableButton(deleteAllBtn);
      actions.push(deleteAllBtn);

      return ui.divH(actions);
    }, true);

    return acc.root;
  }

  // -- Actions --

  private async editCollection(info: MonomerCollectionInfo): Promise<void> {
    const view = MonomerCollectionHandler._view();
    if (view !== null) {
      view.editCollectionPublic(info.name, info.data);
      return;
    }
    const libHelper = await MonomerLibManager.getInstance();
    await libHelper.awaitLoaded();
    await CollectionDialogs.edit(libHelper, libHelper.getMonomerLib(), info, () => this.refreshView());
  }

  private async deleteCollection(info: MonomerCollectionInfo): Promise<void> {
    CollectionDialogs.delete(await MonomerLibManager.getInstance(), [info], () => this.refreshView());
  }

  private async duplicateCollection(info: MonomerCollectionInfo): Promise<void> {
    await CollectionDialogs.duplicate(await MonomerLibManager.getInstance(), info, () => this.refreshView());
  }

  private async mergeCollections(items: MonomerCollectionInfo[]): Promise<void> {
    await CollectionDialogs.merge(await MonomerLibManager.getInstance(), items, () => this.refreshView());
  }

  private async deleteMultipleCollections(items: MonomerCollectionInfo[]): Promise<void> {
    CollectionDialogs.delete(await MonomerLibManager.getInstance(), items, () => this.refreshView());
  }

  private async refreshView(): Promise<void> {
    await MonomerCollectionHandler._view()?.refresh();
  }

  /** The open Monomer Collections view, if any — found from the view's root: the platform hands out a new
   * wrapper for a view each time, so nothing stored on one is there the next. */
  private static _view(): MonomerCollectionsView | null {
    const view = Array.from(grok.shell.views).find((v) => v.name === 'Monomer Collections');
    const root = view?.root.querySelector<HTMLElement>('.monomer-collections-view');
    return root ? (Control.forElement(root) as (Control & MonomerCollectionsView) | undefined) ?? null : null;
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
    const isOwner = CollectionDialogs.canEdit(info.data, currentUser);

    menu.item('Edit', () => { if (isOwner) handler.editCollection(info); }, null,
      {isEnabled: () => CollectionDialogs.notEditable(info.data, currentUser)});
    menu.item('Delete', () => { if (isOwner) handler.deleteCollection(info); }, null,
      {isEnabled: () => CollectionDialogs.notEditable(info.data, currentUser)});
    menu.item('Duplicate', () => handler.duplicateCollection(info));
  }

  /** Builds a context menu for multiple collections. */
  static buildMultiContextMenu(menu: DG.Menu, items: MonomerCollectionInfo[], handler: MonomerCollectionHandler): void {
    const currentUser = DG.User.current().login;
    const allOwned = items.every((c) => CollectionDialogs.canEdit(c.data, currentUser));

    menu.item('Merge', () => handler.mergeCollections(items));
    menu.item(`Delete ${items.length} Collections`, () => { if (allOwned) handler.deleteMultipleCollections(items); },
      null, {isEnabled: () => MonomerCollectionHandler._notAllEditable(items, currentUser)});
  }

  /** Why one of `items` keeps the selection from being changed, naming it, or null. */
  private static _notAllEditable(items: MonomerCollectionInfo[], user: string): string | null {
    const locked = items.find((c) => !CollectionDialogs.canEdit(c.data, user));
    return locked === undefined ? null : `Only ${locked.data.updatedBy} can change ${locked.displayName}`;
  }
}
