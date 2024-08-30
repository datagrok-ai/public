import * as ui from "datagrok-api/ui";
import {Subject} from "rxjs";
import {debounceTime} from "rxjs/operators";
import * as grok from "datagrok-api/grok";
import {DbEntityType} from "./entity";
import {Exp, IFilterDescription, TFilter} from "./cruddy";
import * as DG from "datagrok-api/dg";

export abstract class CruddyFilter {
  root: HTMLDivElement = ui.divV([]);
  entityType: DbEntityType;
  filter: IFilterDescription;
  onChanged: Subject<any> = new Subject<any>();

  abstract getCondition(): TFilter;

  static mergeConditions(conditions: TFilter[]): TFilter {
    const query = {};
    for (const f of conditions)
      Object.assign(query, f);
    return query;
  }

  protected constructor(entityType: DbEntityType, filter: IFilterDescription) {
    this.entityType = entityType;
    this.filter = filter;
  }

  static create(entityType: DbEntityType, filter: IFilterDescription): CruddyFilter {
    switch (filter.type) {
      case 'distinct':
        return new CruddyFilterCategorical(entityType, filter);
      case 'range':
        return new CruddyFilterRange(entityType, filter);
      case 'combo':
        return new CruddyFilterCombo(entityType, filter);
      case 'expression':
        return new CruddyFilterExpression(entityType, filter);
    }

    throw `Unknown filter type: ${filter.type}`;
  }
}

export class CruddyFilterCategorical extends CruddyFilter {
  choices?: DG.InputBase;

  constructor(entityType: DbEntityType, filter: IFilterDescription) {
    super(entityType, filter);
    const tableColumn = filter.column.split('.');
    const tableName = tableColumn.length == 1 ? entityType.table.name : tableColumn[0];

    const sql =
      `select ${filter.column}, count(${filter.column}) 
from ${tableName} 
group by ${filter.column}
order by count(${filter.column}) desc`;

    entityType.crud
      .query(sql)
      //.read({}, { distinct: true, columnNames: [filter.column]})
      .then((df) => {
        this.choices = ui.input.multiChoice('values', {value: [], items: df.columns.byIndex(0).toList()});
        this.choices.onChanged.subscribe(() => this.onChanged.next(this));
        this.root.appendChild(ui.h2(filter.column));
        this.root.appendChild(this.choices.input);

        // counts
        const counts = df.columns.byIndex(1);
        for (let i = 0; i < df.rowCount; i++) {
          const div = this.choices.input.children[i] as HTMLDivElement;
          const percent = Math.round(counts.getNumber(i) * 100 / counts.stats.max);
          div.style.background = `linear-gradient(to right, #eff7fa ${percent}%, white ${percent}%)`;
          div.appendChild(ui.divText(`${counts.getNumber(i)}`, 'cruddy-category-count'));
        }
      });
  }

  getCondition(): TFilter {
    if (!this.choices) return {};
    const choices = (this.choices!.value as Array<any>);
    if (choices.length == 0) return {};

    return {
      [`${this.filter.column} in (${choices.map((x) => `'${x}'`).join(', ')})`]: null
    };
  }
}

export class CruddyFilterCombo extends CruddyFilter {
  choices?: DG.ChoiceInput<string>;

  constructor(entityType: DbEntityType, filter: IFilterDescription) {
    super(entityType, filter);
    entityType.crud
      .read({}, {distinct: true, columnNames: [filter.column]})
      .then((df) => {
        const items = df.columns.byIndex(0).toList();
        this.choices = ui.input.choice('values', {value: null, items: items, nullable: true});
        this.choices.onChanged.subscribe(() => this.onChanged.next(this));
        this.root.appendChild(ui.h2(filter.column));
        this.root.appendChild(this.choices.input);
      });
  }

  getCondition(): TFilter {
    return this.choices?.value == null ? {} : {[this.filter.column]: this.choices.value};
  }
}

export class CruddyFilterExpression extends CruddyFilter {

  operation = ui.input.choice('operator', {value: '', items: ['']});
  colInput = ui.input.choice('column', {value: this.entityType.columns[0].name, items: this.entityType.columns.map((c) => c.name)});
  input = ui.input.string('input', {value: ''});

  constructor(entityType: DbEntityType, filter: IFilterDescription) {
    super(entityType, filter);
    this.root.className = 'd4-flex-row';
    this.root.appendChild(this.colInput.input);
    this.root.appendChild(this.operation.input);
    this.root.appendChild(this.input.input);

    const refresh = () => {
      this.operation.items = Exp.typeOperators[this.entityType.getColumn(this.colInput.value!).type];
    };

    this.colInput.onChanged.subscribe(() => {
      refresh();
      this.onChanged.next(this);
    });
    this.operation.onChanged.subscribe(() => this.onChanged.next(this));
    this.input.onChanged.subscribe(() => this.onChanged.next(this));

    refresh();
  }

  getCondition(): TFilter {
    return Exp.getCondition(this.colInput.value!, this.operation.value!, this.input.value);
  }
}

export class CruddyFilterRange extends CruddyFilter {
  slider = ui.rangeSlider(0, 1, 0, 1, false, 'thin_barbell');

  constructor(entityType: DbEntityType, filter: IFilterDescription) {
    super(entityType, filter);
    this.slider.root.querySelector('svg')!.style.height = '20px';
    const minCol = `min(${filter.column}) as min`;
    const maxCol = `max(${filter.column}) as max`;
    entityType.crud
      .read({}, {columnNames: [minCol, maxCol]})
      .then((df) => {
        this.slider.setValues(df.get('min', 0), df.get('max', 0), df.get('min', 0), df.get('max', 0));
        this.slider.onValuesChanged.subscribe((_) => this.onChanged.next(this));
        this.root.appendChild(ui.h2(filter.column));
        this.root.appendChild(this.slider.root);
      });
  }

  getCondition(): TFilter {
    return {
      [`(${this.filter.column} >= ${this.slider.min} and ${this.filter.column} <= ${this.slider.max})`]: null
    };
  }
}

export class CruddyFilterHost {
  root: HTMLDivElement = ui.divV([]);
  onChanged: Subject<TFilter> = new Subject<TFilter>();
  filters: CruddyFilter[] = [];

  constructor() {
    this.onChanged
      .pipe(debounceTime(100))
      .subscribe((_) => grok.shell.info('Filter changed: ' + JSON.stringify(this.getCondition())));
  }

  getCondition(): TFilter {
    return CruddyFilter.mergeConditions(this.filters.map((f) => f.getCondition()));
  }

  init(entityType: DbEntityType) {
    for (const f of entityType.filters) {
      const filter = CruddyFilter.create(entityType, f);
      filter.onChanged.subscribe((x) => this.onChanged.next(x));
      this.filters.push(filter);
      this.root.appendChild(filter.root);
    }
  }
}