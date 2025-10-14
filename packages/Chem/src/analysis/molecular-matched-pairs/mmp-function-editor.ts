import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ColumnInputOptions} from '@datagrok-libraries/utils/src/type-declarations';
import {SCALING_METHODS} from './mmp-viewer/mmp-constants';

export enum MmpDiffTypes {
  delta = 'delta',
  ratio = 'ratio'
}

export type MMPFunctionParams = {
  table: DG.DataFrame,
  molecules: DG.Column,
  activities: DG.ColumnList,
  diffTypes: MmpDiffTypes[],
  scalings: SCALING_METHODS[],
  fragmentCutoff: number
}

export class MmmpFunctionEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  colInput!: DG.InputBase<DG.Column | null>;
  colInputRoot!: HTMLElement;
  activitiesInput!: DG.InputBase;
  activitiesInputRoot = ui.div();
  cutoffInput = ui.input.float('Cutoff', {value: 0.4});
  activitiesParams: {[key: string]: {deltaType: string, scaling: string}} = {};
  activitiesParamsDiv = ui.divV([]);

  constructor() {
    this.tableInput =
        ui.input.table('Table', {value: grok.shell.tv.dataFrame, items: grok.shell.tables, onValueChanged: () => {
          this.onTableInputChanged();
        }});
    this.onTableInputChanged();
  }

  onTableInputChanged() {
    this.activitiesParams = {};
    ui.empty(this.activitiesInputRoot);
    const numericalColumns = Array.from(this.tableInput.value!.columns.numerical);
    this.activitiesInput = ui.input.columns('Activities', {
      table: this.tableInput.value!,
      filter: (col: DG.Column) => numericalColumns.includes(col),
      onValueChanged: () => {
        ui.empty(this.activitiesParamsDiv);
        const df = this.tableInput.value!;
        const colNames: string[] = this.activitiesInput.value.map((it: DG.Column) => it.name);
        for (const key of (Object.keys(this.activitiesParams))) {
          if (!colNames.includes(key))
            delete this.activitiesParams[key];
        }
        colNames.forEach((it) => {
          const typeChoice = ui.input.choice(it, {
            items: Object.values(MmpDiffTypes),
            value: this.activitiesParams[it]?.deltaType ?? MmpDiffTypes.delta,
            onValueChanged: () => {
              this.activitiesParams[it].deltaType = typeChoice.value!;
            },
          });
          const scalingChoice = ui.input.choice('Scaling', {
            items: df.col(it)!.stats.min < 0 ? [SCALING_METHODS.NONE] :
              [SCALING_METHODS.NONE, SCALING_METHODS.LG, SCALING_METHODS.MINUS_LG],
            value: this.activitiesParams[it]?.scaling ?? SCALING_METHODS.NONE,
            onValueChanged: () => {
              this.activitiesParams[it].scaling = scalingChoice.value!;
            },
          });
          this.activitiesParams[it as string] = {deltaType: typeChoice.value!, scaling: scalingChoice.value!};
          this.activitiesParamsDiv.append(ui.divH([typeChoice.root, scalingChoice.root]));
        });
      },
    } as ColumnInputOptions);
    this.activitiesInputRoot.append(this.activitiesInput.root);
    this.regenerateColInput();
  }

  private getColInput() {
    const firstSupportedColumn = this.tableInput.value?.columns.toList()
      .find((col) => col.semType === DG.SEMTYPE.MOLECULE) ?? undefined;
    const input = ui.input.column('Column', {
      table: this.tableInput.value!,
      filter: (col) => col.semType === DG.SEMTYPE.MOLECULE,
      value: firstSupportedColumn});
    if (!this.colInputRoot)
      this.colInputRoot = input.root;
    return input;
  }

  private regenerateColInput() {
    let flag = false;
    if (this.colInputRoot) {
      flag = true;
      ui.empty(this.colInputRoot);
    }
    this.colInput = this.getColInput();
    if (flag)
      Array.from(this.colInput.root.children).forEach((child) => this.colInputRoot.append(child));
  }

  public getEditor(): HTMLElement {
    return ui.div([
      this.tableInput,
      this.colInputRoot,
      this.activitiesInputRoot,
      this.activitiesParamsDiv,
      this.cutoffInput,
    ], {style: {minWidth: '320px'}});
  }

  getDiffTypes(): MmpDiffTypes[] {
    const arr: MmpDiffTypes[] = [];
    this.activitiesInput!.value.forEach((col: DG.Column) => {
      arr.push(this.activitiesParams[col.name].deltaType as MmpDiffTypes);
    });
    return arr;
  }

  getScalings(): SCALING_METHODS[] {
    const arr: SCALING_METHODS[] = [];
    this.activitiesInput!.value.forEach((col: DG.Column) => {
      arr.push(this.activitiesParams[col.name].scaling as SCALING_METHODS);
    });
    return arr;
  }

  public getParams(): MMPFunctionParams {
    return {
      table: this.tableInput.value!,
      molecules: this.colInput.value!,
      activities: this.activitiesInput!.value,
      diffTypes: this.getDiffTypes(),
      scalings: this.getScalings(),
      fragmentCutoff: this.cutoffInput.value!,
    };
  }
}
