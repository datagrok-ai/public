import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Subject} from 'rxjs';

export type ItemType = {[_:string]: any};
export type ItemsGridOptions = {
    removeButtonTooltip?: string,
    addButtonTooltip?: string,
    allowRemove?: boolean,
    allowAdd?: boolean,
    horizontalInputNames?: boolean,
    newItemFunction?: () => ItemType,
};

/** Editor for predefined properties in a grid form. Supports adding, removing and editing of props. */
export class ItemsGrid {
  private _root: HTMLDivElement;
  private _items: {[key: string]: any}[];

  /** Observable for when the added item is change */
  public onItemChanged: Subject<{item:ItemType, fieldName: string}> = new Subject();
  /** Observable for when the adding item is change */
  public onAddingItemChanged: Subject<{item:ItemType, fieldName: string}> = new Subject();
  /** Observable for when the item is added */
  public onItemAdded: Subject<ItemType> = new Subject();
  /** Observable for when the item is removed */
  public onItemRemoved: Subject<ItemType> = new Subject();

  properties: DG.Property[];
  options: ItemsGridOptions =
    {removeButtonTooltip: 'Remove item', addButtonTooltip: 'Add item', allowAdd: true, allowRemove: true};
    /**
     * Creates an instance of ItemsGrid which is a grid form for editing properties.
     * @param {DG.Property[]}properties - properties to be edited
     * @param {ItemTypep[]}items - items to be edited
     * @param {ItemsGridOptions}options - options for the grid
     * Options include:
     * removeButtonTooltip - tooltip for the remove button - default is 'Remove item'
     * addButtonTooltip - tooltip for the add button - default is 'Add item'
     * allowRemove - whether to allow removing items - default is true
     * allowAdd - whether to allow adding items - default is true
     * horizontalInputNames - Show the input names next to intput - default is false, names are shown above the grid.
     * newItemFunction - function that returns new item - default is undefind, in which case an empty object is created.
     */
  constructor(properties: DG.Property[], items: ItemType[] = [], options: ItemsGridOptions = {}) {
    this.properties = properties;
    this.options = {...this.options, ...options};
    this._root = ui.divV([],
      {style: {
        display: 'grid', gridTemplateColumns: `repeat(${this.properties.length}, 1fr)`,
        alignItems: 'center', gap: '12px'}});
    this._items = items;
    this.render();
  }

  get root(): HTMLDivElement {
    return this._root;
  }

  get items(): {[key: string]: any}[] {
    return this._items;
  }

  private render(): void {
    ui.empty(this._root);
    if (!this.options.horizontalInputNames) {
      for (const prop of this.properties) {
        const header = ui.divText(prop.name);
        header.style.fontWeight = 'bold';
        this._root.appendChild(header);
      }
    }
    for (const item of this._items) {
      const editors = this.getItemDiv(item);
      // if (editors.length !== this.properties.length + 1)
      //   editors.push(ui.div());
      for (const editor of editors)
        this._root.appendChild(editor);
    }

    if (this.options.allowAdd) {
      const newEditor = this.getItemDiv(undefined, true);
      for (const editor of newEditor)
        this._root.appendChild(editor);
    }
  }

  private getItemDiv(item: ItemType = {}, isAdding?: boolean): HTMLElement[] {
    const editors: HTMLElement[] = [];

    const inputsMap: {[_: string]: DG.InputBase} = {};
    let lastInput: DG.InputBase | null = null;
    for (const prop of this.properties) {
      if (item[prop.name] === undefined)
        item[prop.name] = null; // needed for date editor, it can not handle undefined
      const input = ui.input.forProperty(prop, item);
      editors.push(this.options.horizontalInputNames ? input.root : this.hideLabel(input.root));
      if (prop.propertyType !== DG.TYPE.BOOL)
        input.input.style.width = '100%';
      inputsMap[prop.name] = input;
      input.onChanged(() => {
        isAdding ? this.onAddingItemChanged.next({item, fieldName: prop.name}) :
          this.onItemChanged.next({item, fieldName: prop.name});
      });
      input.root.style.alignItems = 'center'; // needed for molecule editor
      lastInput = input;
    }

    let companionButton: HTMLElement | null = null;
    if (isAdding) {
      const addButton = ui.icons.add(() => {
        const newItem: ItemType = this.options.newItemFunction ? this.options.newItemFunction() : {};
        Object.keys(inputsMap).forEach((pName) => {
          newItem[pName] = inputsMap[pName].value;
        });
        this.onItemAdded.next(newItem);
        this._items.push(newItem);
        this.render();
      }, this.options.addButtonTooltip);
      companionButton = addButton;
    } else {
      if (!this.options.allowRemove) return editors;
      const removeButton = ui.icons.delete(() => {
        this.onItemRemoved.next(item);
        this._items.splice(this._items.indexOf(item), 1);
        this.render();
      }, this.options.removeButtonTooltip);
      companionButton = removeButton;
    }

    //editors.push(companionButton);
    lastInput && lastInput.addOptions(companionButton);
    companionButton.style.color = '#2083d5';
    return editors;
  }
  //needed for color input
  private hideLabel(el: HTMLElement): HTMLElement {
    el.removeChild(el.getElementsByTagName('label')[0]);
    return el;
  }
}

/** Usage example:
 * const props = [DG.Property.js('name2', DG.TYPE.STRING),
      DG.Property.js('name3', DG.TYPE.NUM),
      DG.Property.fromOptions({name: 'name4', type: DG.TYPE.DATE_TIME, nullable: true})];
    const itemsGrid = new ItemsGrid(props);
 * the component has a root property that can be added to the page
    ui.dialog('dialog').add(itemsGrid.root).show();
    at any point, filled fields can be accessed via itemsGrid.items
 */
