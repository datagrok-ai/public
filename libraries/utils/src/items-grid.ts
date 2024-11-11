import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Subject, Observable} from 'rxjs';

export type ItemType = {[_:string]: any};
export type ItemsGridOptions = {
    removeButtonTooltip?: string,
    addButtonTooltip?: string,
    allowRemove?: boolean,
    allowAdd?: boolean,
    horizontalInputNames?: boolean,
    newItemFunction?: () => ItemType,
    customInputs?: {[key: string]: (item: ItemType) => InputType},
    validators?: {[key: string]: (value: string) => string | null},
    customLabels?: {[key: string]: string},
    inputOptions?: {[key: string]: HTMLElement}
};

export type InputType = {
    input?: HTMLElement,
    root: HTMLElement,
    value: any,
    onChanged: Observable<any>,
    addOptions?: (el: HTMLElement) => void,
    addValidator: (validator: (value: string) => string | null) => void
}

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
        alignItems: 'center', gap: '12px'}, classes: 'ui-items-grid'});
    this._items = items;
    this.render();
  }

  get root(): HTMLDivElement {
    return this._root;
  }

  get items(): {[key: string]: any}[] {
    return this._items;
  }

  set items(value: {[key: string]: any}[]) {
    this._items = value;
    this.render();
  }

  public addItem(item: ItemType, notify = true): void {
    this._items.push(item);
    this.render();
    if (notify)
      this.onItemAdded.next(item);
  }

  public removeItem(item: ItemType, notify = true): void {
    if (this._items.indexOf(item) === -1)
      return;
    this._items.splice(this._items.indexOf(item), 1);
    this.render();
    if (notify)
      this.onItemRemoved.next(item);
  }

  public removeAtIndex(index: number, notify = true): void {
    if (index < 0 || index >= this._items.length)
      return;
    const removed = this._items.splice(index, 1);
    this.render();
    if (notify)
      this.onItemRemoved.next(removed[0]);
  }

  public removeAllItems(): void {
    this._items = [];
    this.render();
  }

  public render(): void {
    ui.empty(this._root);
    if (!this.options.horizontalInputNames) {
      for (const prop of this.properties) {
        let label = prop.caption ?? prop.name;
        if (this.options.customLabels?.[prop.name])
          label = this.options.customLabels[prop.name];
        const header = ui.divText(label);
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
      const newEditor = this.getItemDiv(this.options.newItemFunction?.() ?? undefined, true);
      for (const editor of newEditor)
        this._root.appendChild(editor);
    }
  }

  private addingItemInputs: {[key: string]: InputType} = {};

  public get addingItem(): ItemType {
    const res: ItemType = {};
    Object.keys(this.addingItemInputs).forEach((pName) => {
      res[pName] = this.addingItemInputs[pName].value;
    });
    return res;
  }

  private getItemDiv(item: ItemType = {}, isAdding?: boolean): HTMLElement[] {
    const editors: HTMLElement[] = [];

    const inputsMap: {[_: string]: InputType} = {};
    let lastInput: InputType | null = null;
    for (const prop of this.properties) {
      if (item[prop.name] === undefined)
        item[prop.name] = null; // needed for date editor, it can not handle undefined
      const input = this.options.customInputs?.[prop.name] ? this.options.customInputs[prop.name](item) :
        ui.input.forProperty(prop, item);
      const validator = this.options.validators?.[prop.name];
      if (validator)
        input.addValidator(validator);

      editors.push(this.options.horizontalInputNames ? input.root : this.hideLabel(input.root));
      if (prop.propertyType !== DG.TYPE.BOOL && prop.name.toLowerCase() !== 'color')
        input.input && (input.input.style.width = '100%');
      inputsMap[prop.name] = input;
      input.onChanged.subscribe(() => {
        item[prop.name] = input.value;
        isAdding ? this.onAddingItemChanged.next({item, fieldName: prop.name}) :
          this.onItemChanged.next({item, fieldName: prop.name});
      });
      input.root.style.alignItems = 'center'; // needed for molecule editor
      lastInput = input;
    }

    let companionButton: HTMLElement | null = null;
    if (isAdding) {
      this.addingItemInputs = inputsMap;
      const addButton = ui.icons.add(() => {
        const newItem: ItemType = this.options.newItemFunction ? this.options.newItemFunction() : {};
        Object.keys(inputsMap).forEach((pName) => {
          newItem[pName] = inputsMap[pName].value;
        });
        this._items.push(newItem);
        this.onItemAdded.next(newItem);
        this.render();
      }, this.options.addButtonTooltip);
      companionButton = addButton;
    } else {
      if (!this.options.allowRemove) return editors;
      const removeButton = ui.icons.delete(() => {
        this._items.splice(this._items.indexOf(item), 1);
        this.onItemRemoved.next(item);
        this.render();
      }, this.options.removeButtonTooltip);
      companionButton = removeButton;
    }

    //editors.push(companionButton);
    lastInput && lastInput.addOptions ? lastInput.addOptions(companionButton) :
      lastInput?.root.appendChild(companionButton);
    companionButton.style.color = '#2083d5';
    return editors;
  }
  //needed for color input
  private hideLabel(el: HTMLElement): HTMLElement {
    el.getElementsByTagName('label')[0] && el.removeChild(el.getElementsByTagName('label')[0]);
    return el;
  }

  hasErrors(): boolean {
    return this._root.querySelectorAll('.d4-invalid').length > 0;
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
