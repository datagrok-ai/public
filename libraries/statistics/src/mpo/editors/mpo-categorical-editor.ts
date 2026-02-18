import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';
import {CategoricalDesirability} from '../mpo';

export class MpoCategoricalEditor {
  root: HTMLDivElement;
  onChanged = new Subject<CategoricalDesirability>();
  supportsModeDialog: boolean = true;

  private _prop: CategoricalDesirability;
  private design: boolean;
  private showControls: boolean;
  private desirabilityInputs: DG.InputBase[] = [];
  private form: HTMLElement | null = null;
  private columnCategories: string[] | null = null;

  constructor(prop: CategoricalDesirability, design = false, showControls = false) {
    this._prop = prop;
    this.design = design;
    this.showControls = showControls;
    this.root = ui.divV([], 'statistics-mpo-cat-editor');
    this.buildForm();
  }

  private buildForm(): void {
    const categories = this._prop.categories ?? [];

    if (this.form) {
      ui.empty(this.root);
      this.form = null;
      this.desirabilityInputs = [];
    }

    if (this.design)
      this.buildDesignForm(categories);
    else
      this.buildViewForm(categories);
  }

  private createDesirabilityInput(
    label: string,
    cat: {name: string; desirability: number},
  ): DG.InputBase {
    const input = ui.input.float(label, {value: cat.desirability ?? 0.5, min: 0, max: 1, format: '#0.000',
      onValueChanged: (value: number) => {
        cat.desirability = value;
        this.onChanged.next(this._prop);
      },
    });
    this.desirabilityInputs.push(input);
    return input;
  }

  private buildViewForm(categories: {name: string; desirability: number}[]): void {
    for (const cat of categories)
      this.createDesirabilityInput(cat.name, cat);

    this.form = ui.form(this.desirabilityInputs);
    this.root.append(this.form);
  }

  private buildDesignForm(categories: {name: string; desirability: number}[]): void {
    const rows = categories.map((cat, idx) => this.buildCategoryRow(cat, idx));

    this.form = ui.divV(rows);
    this.root.append(this.form);
  }

  private buildCategoryRow(
    cat: {name: string; desirability: number},
    idx: number,
  ): HTMLElement {
    const nameInput = this.columnCategories ?
      ui.input.choice('', {items: this.columnCategories, nullable: true, value: cat.name || null,
        onValueChanged: (v) => {
          cat.name = v ?? '';
          this.onChanged.next(this._prop);
        },
      }) :
      ui.input.string('', {value: cat.name, onValueChanged: (v) => {
        cat.name = v;
        this.onChanged.next(this._prop);
      }});
    nameInput.root.classList.add('statistics-mpo-cat-name');

    const desInput = this.createDesirabilityInput('', cat);
    desInput.root.classList.add('statistics-mpo-cat-desirability');

    const elements: HTMLElement[] = [nameInput.root, desInput.root];

    if (this.showControls) {
      const add = ui.icons.add(() => {
        const newCat = {name: this.columnCategories ? '' : `Category ${this._prop.categories.length + 1}`, desirability: 1};
        this._prop.categories.splice(idx + 1, 0, newCat);
        this.buildForm();
        this.onChanged.next(this._prop);
      });

      const del = ui.icons.delete(() => {
        this._prop.categories.splice(idx, 1);
        this.buildForm();
        this.onChanged.next(this._prop);
      });

      elements.push(ui.divH([add, del], 'statistics-mpo-control-buttons'));
    }

    return ui.divH(elements, 'statistics-mpo-cat-row');
  }

  redrawAll(notify: boolean = true): void {
    this.buildForm();
    if (notify)
      this.onChanged.next(this._prop);
  }

  setChoices(choices: string[]): void {
    this.columnCategories = choices;
    this.buildForm();
  }

  setColumn(col: DG.Column | null): void {
    if (!col)
      return;

    this.columnCategories = col.isCategorical ? [...col.categories] : null;
    const existing = new Map(this._prop.categories?.map((c) => [c.name, c.desirability]) ?? []);
    this._prop.categories = col.categories.map((c: string) => ({
      name: c,
      desirability: existing.get(c) ?? 1,
    }));
    this.buildForm();
  }
}
