import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Subject} from "rxjs";

export class SchemaEditor extends DG.Widget {
  table: PropertyTable;
  allowedTypes: string[] = [DG.TYPE.STRING, DG.TYPE.INT, DG.TYPE.FLOAT, DG.TYPE.BOOL, DG.TYPE.DATE_TIME];

  constructor(options: {properties: DG.PropertyOptions[]}) {
    super(ui.div([]));

    this.table = new PropertyTable({
      items: options.properties,
      mainProperties: [DG.Property.propertyOptions['name']!, DG.Property.propertyOptions['type']!],
      allowAdd: true,
      allowRemove: true
    });

    this.root.appendChild(this.table.root);
  }
}

export type TableFromPropertiesOptions<T = any> = {
  items: T[],
  mainProperties: DG.PropertyOptions[],
  extraProperties?: DG.PropertyOptions[],
  allowAdd?: boolean,
  allowRemove?: boolean,
  createNew?: () => T
}

export class PropertyTable<T = any> extends DG.Widget {
  table?: HTMLTableElement;
  onItemsChanged: Subject<null> = new Subject<null>();
  onItemPropertyChanged: Subject<null> = new Subject<null>();
  onItemAdded: Subject<null> = new Subject<null>();
  onItemRemoved: Subject<null> = new Subject<null>();

  constructor(public data: TableFromPropertiesOptions) {
    super(ui.div([]));
    this.refresh();
  }

  refresh() {
    const createRow = (item: T, i: number) => {
      const elements = this.data.mainProperties
        .map((p) => DG.InputBase.forProperty(DG.Property.fromOptions(p), item).input);
      if (this.data.allowAdd ?? false)
        elements.push(ui.iconFAB('plus', () => {
          const newItem = this.data.createNew ? this.data.createNew() : {};
          this.data.items.splice(this.data.items.indexOf(item) + 1, 0, newItem);
          this.refresh();
          this.onItemAdded.next();
        }));
      if (this.data.allowRemove ?? false)
        elements.push(ui.iconFAB('times', () => {
          this.data.items.splice(this.data.items.indexOf(item), 1);
          this.refresh();
          this.onItemRemoved.next();
        }));
      return elements;
    }

    this.table = ui.table(this.data.items, createRow);
    ui.empty(this.root);
    this.root.appendChild(this.table);
  }
}