import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject} from 'rxjs';
import {PropertyDesirability} from './mpo';

export class MpoCategoricalEditor {
  root = ui.divV([]);
  onChanged = new Subject<PropertyDesirability>();

  private _prop: PropertyDesirability;
  private form = ui.form([]);

  constructor(prop: PropertyDesirability) {
    this._prop = prop;
    // this.root.appendChild(this.form);

    this.buildForm();
  }

  private buildForm() {
    const categories = this._prop.categories ?? [];
    const inputs: DG.InputBase[] = [];

    categories.forEach((cat, idx) => {
      const slider = ui.input.slider(cat.name, {
        value: cat.weight ?? 0.5,
        min: 0,
        max: 1,
        step: 0.01,
        onValueChanged: (v: number) => {
          this._prop.categories![idx].weight = v;
          this.onChanged.next(this._prop);
        },
      });
      inputs.push(slider);

      // const row = ui.divH([
      //   ui.label(cat.name, {style: {width: '70px'}}),
      //   slider.root,
      // ], {style: {alignItems: 'center', gap: '10px'}});

      // this.form.appendChild(row);
    });
    this.form = ui.form(inputs);
    this.root.append(this.form);
  }

  drawBars(values?: number[]) {};
}

