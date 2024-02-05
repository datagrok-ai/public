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
    entityType.crud
      .query(`select ${filter.column}, count(${filter.column}) from ${entityType.table.name} group by ${filter.column}`)
      //.read({}, { distinct: true, columnNames: [filter.column]})
      .then((df) => {
        this.choices = ui.multiChoiceInput('values', [], df.columns.byIndex(0).toList());
        this.choices.onChanged(() => this.onChanged.next(this));
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
        this.choices = ui.choiceInput('values', null, items, null, {nullable: true});
        this.choices.onChanged(() => this.onChanged.next(this));
        this.root.appendChild(ui.h2(filter.column));
        this.root.appendChild(this.choices.input);
      });
  }

  getCondition(): TFilter {
    return this.choices?.value == null ? {} : {[this.filter.column]: this.choices.value};
  }
}

export class CruddyFilterExpression extends CruddyFilter {

  operation = ui.choiceInput('operator', '', ['']);
  colInput = ui.choiceInput('column', this.entityType.columns[0].name, this.entityType.columns.map((c) => c.name));
  input = ui.stringInput('input', '');

  constructor(entityType: DbEntityType, filter: IFilterDescription) {
    super(entityType, filter);
    this.root.className = 'd4-flex-row';
    this.root.appendChild(this.colInput.input);
    this.root.appendChild(this.operation.input);
    this.root.appendChild(this.input.input);

    const refresh = () => {
      this.operation.items = Exp.typeOperators[this.entityType.getColumn(this.colInput.value!).type];
    };

    this.colInput.onChanged(() => {
      refresh();
      this.onChanged.next(this);
    });
    this.operation.onChanged(() => this.onChanged.next(this));
    this.input.onChanged(() => this.onChanged.next(this));

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