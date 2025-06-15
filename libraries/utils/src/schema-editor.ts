import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {merge, Observable, Subject} from "rxjs";
import {IProperty, ValueMatcher} from "datagrok-api/dg";

export class SchemaEditor extends DG.Widget {
  table: PropertyTable;
  allowedTypes: string[] = [DG.TYPE.STRING, DG.TYPE.INT, DG.TYPE.FLOAT, DG.TYPE.BOOL, DG.TYPE.DATE_TIME];
  properties: DG.IProperty[] = [];

  constructor(options: {properties: DG.IProperty[]}) {
    super(ui.div([]));
    this.properties = options.properties;

    const typeProp: IProperty = {
      ...DG.Property.propertyOptions['type'],
      choices: this.allowedTypes
    };

    this.table = new PropertyTable({
      items: options.properties,
      mainProperties: [DG.Property.propertyOptions['name']!, typeProp],
      allowAdd: true,
      allowRemove: true
    });

    this.root.appendChild(this.table.root);
  }
}

export type TableFromPropertiesOptions<T = any> = {
  items: T[],
  mainProperties: DG.IProperty[],
  extraProperties?: DG.IProperty[],
  allowAdd?: boolean,
  allowRemove?: boolean,
  createNew?: () => T
}

export class PropertyTable<T = any> extends DG.Widget {
  table?: HTMLTableElement;

  onSelected: Subject<T> = new Subject();
  onItemChanged: Subject<T> = new Subject();
  onItemAdded: Subject<T> = new Subject();
  onItemRemoved: Subject<T> = new Subject();

  /** Any change (ItemChanged, ItemAdded, ItemRemoved) */
  get onChanged(): Observable<T> {
    return merge(this.onItemChanged, this.onItemAdded, this.onItemRemoved);
  }

  constructor(public data: TableFromPropertiesOptions) {
    super(ui.div([]));
    this.refresh();
  }

  refresh() {
    const createRow = (item: T, i: number) => {
      const elements = this.data.mainProperties
        .map((p) => {
          const input = DG.InputBase.forProperty(DG.Property.fromOptions(p), item);
          for (const v of p.validators ?? [])
            input.addValidator(ValueMatcher.forType(p.type!, v).validate);
          input.onChanged.subscribe(() => this.onItemChanged.next(item));
          input.input.onclick = () => this.onSelected.next(item);
          return input.input;
        });

      if (this.data.allowAdd ?? false)
        elements.push(ui.iconFAB('plus', () => {
          const newItem = this.data.createNew ? this.data.createNew() : {};
          this.data.items.splice(this.data.items.indexOf(item) + 1, 0, newItem);
          this.refresh();
          this.onItemAdded.next(newItem);
        }));
      if (this.data.allowRemove ?? false)
        elements.push(ui.iconFAB('times', () => {
          this.data.items.splice(this.data.items.indexOf(item), 1);
          this.refresh();
          this.onItemRemoved.next(item);
        }));
      return elements;
    }

    this.table = ui.table(this.data.items, createRow);
    ui.empty(this.root);
    this.root.appendChild(this.table);
  }
}


export class PropertyValidator {
  static NOT_EMPTY = 'not empty';

  static validators: {[name: string]: (value: any) => string | null} = {
    NOT_EMPTY: (s: any) => s == null || s === '' ? "Can't be empty" : null
  }

  static getValidator(obj: any, prop: DG.IProperty, expression: string): ((value: any) => string | null) {
    if (this.validators[expression] !== null)
      return this.validators[expression];

    if (ValueMatcher.supportedTypes.includes(prop.type!))
      return ValueMatcher.forType(prop.type!, expression).validate;

    return (_) => null;
  }

  static registerValidator(name: string, validator: ((x: any) => string | null)) {
    this.validators[name] = validator;
  }
}