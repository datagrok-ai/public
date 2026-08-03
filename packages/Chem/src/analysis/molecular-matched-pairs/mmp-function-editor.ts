import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, Subject} from 'rxjs';
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
  fragmentCutoff: number,
  runOnFilteredData: boolean
}

/** FuncCall editor for `Chem:mmpAnalysis`: edits the live funccall inputs (table, molecules,
 * activities, per-activity diffTypes/scalings, fragmentCutoff, runOnFilteredData); the platform
 * hosts it in a dialog and runs the call on OK. */
export class MmmpFunctionEditor extends DG.FuncCallEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  colInput!: DG.InputBase<DG.Column | null>;
  colInputRoot!: HTMLElement;
  activitiesInput!: DG.InputBase;
  activitiesInputRoot = ui.div();
  cutoffInput = ui.input.float('Cutoff', {value: 0.4, min: 0, max: 1, nullable: false});
  activitiesParams: {[key: string]: {deltaType: string, scaling: string}} = {};
  activitiesParamsDiv = ui.divV([]);
  runOnFilteredDataInput: DG.InputBase<boolean>;
  private inputChangedSubject: Subject<any> = new Subject<any>();

  constructor(private funcCall: DG.FuncCall) {
    const root = ui.div([]);
    super(root);
    this.tableInput =
        ui.input.table('Table', {value: funcCall.inputs['table'] ?? grok.shell.tv.dataFrame,
          items: grok.shell.tables, onValueChanged: () => {
            this.onTableInputChanged();
            this.syncCall();
          }});
    this.onTableInputChanged();
    this.runOnFilteredDataInput = ui.input.bool('Run On Filtered Data', {value: true,
      tooltipText: 'If checked, MMP analysis will run on filtered data', onValueChanged: () => this.syncCall()});
    this.cutoffInput.onChanged.subscribe(() => this.syncCall());
    root.append(this.getEditor());
    this.syncCall();
  }

  /** Writes the current UI selection into the live funccall and notifies the platform. */
  private syncCall(): void {
    this.funcCall.inputs['table'] = this.tableInput.value;
    this.funcCall.inputs['molecules'] = this.colInput?.value;
    this.funcCall.inputs['activities'] = this.activitiesInput?.value;
    this.funcCall.inputs['diffTypes'] = this.getDiffTypes();
    this.funcCall.inputs['scalings'] = this.getScalings();
    this.funcCall.inputs['fragmentCutoff'] = this.cutoffInput.value;
    this.funcCall.inputs['runOnFilteredData'] = this.runOnFilteredDataInput?.value ?? true;
    this.inputChangedSubject.next(null);
  }

  onTableInputChanged() {
    this.activitiesParams = {};
    ui.empty(this.activitiesInputRoot);
    const numericalColumns = Array.from(this.tableInput.value!.columns.numerical);
    this.activitiesInput = ui.input.columns('Activities', {
      table: this.tableInput.value!,
      filter: (col: DG.Column) => numericalColumns.includes(col),
      nullable: false,
      onValueChanged: () => {
        this.updateActivitiesParamsDiv();
        this.syncCall();
      },
    } as ColumnInputOptions);
    this.activitiesInputRoot.append(this.activitiesInput.root);
    this.regenerateColInput();
  }

  /** Rebuilds the per-activity delta/scaling choice rows from the current activities selection,
   * preserving already-chosen mappings in {@link activitiesParams}. Called on user selection and
   * explicitly on history restore (a programmatic value set does not fire onValueChanged). */
  private updateActivitiesParamsDiv(): void {
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
          this.syncCall();
        },
      });
      const scalingChoice = ui.input.choice('Scaling', {
        items: df.col(it)!.stats.min < 0 ? [SCALING_METHODS.NONE] :
          [SCALING_METHODS.NONE, SCALING_METHODS.LG, SCALING_METHODS.MINUS_LG],
        value: this.activitiesParams[it]?.scaling ?? SCALING_METHODS.NONE,
        onValueChanged: () => {
          this.activitiesParams[it].scaling = scalingChoice.value!;
          this.syncCall();
        },
      });
      this.activitiesParams[it as string] = {deltaType: typeChoice.value!, scaling: scalingChoice.value!};
      this.activitiesParamsDiv.append(ui.divH([typeChoice.root, scalingChoice.root]));
    });
  }

  private getColInput() {
    const firstSupportedColumn = this.tableInput.value?.columns.toList()
      .find((col) => col.semType === DG.SEMTYPE.MOLECULE) ?? undefined;
    const input = ui.input.column('Column', {
      table: this.tableInput.value!,
      filter: (col) => col.semType === DG.SEMTYPE.MOLECULE,
      value: firstSupportedColumn,
      nullable: false,
      onValueChanged: () => this.syncCall(),
    });
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
      this.runOnFilteredDataInput,
    ], {style: {minWidth: '320px'}});
  }

  getDiffTypes(): MmpDiffTypes[] {
    const arr: MmpDiffTypes[] = [];
    this.activitiesInput?.value?.forEach((col: DG.Column) => {
      arr.push(this.activitiesParams[col.name].deltaType as MmpDiffTypes);
    });
    return arr;
  }

  getScalings(): SCALING_METHODS[] {
    const arr: SCALING_METHODS[] = [];
    this.activitiesInput?.value?.forEach((col: DG.Column) => {
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
      runOnFilteredData: this.runOnFilteredDataInput.value!,
    };
  }

  get isValid(): boolean {
    // Pure read only: the platform re-evaluates isValid on every onInputChanged emission, so
    // emitting from here (via syncCall) would spin an infinite loop. Inputs are kept current by
    // the onValueChanged handlers, which is the last thing to fire before the platform reads this.
    return this.tableInput.value != null && this.colInput?.value != null &&
      (this.activitiesInput?.value?.length ?? 0) >= 1 &&
      this.cutoffInput.value != null && this.cutoffInput.value >= 0 && this.cutoffInput.value <= 1;
  }

  getHistoryString(): string {
    return JSON.stringify({
      activities: (this.activitiesInput?.value ?? []).map((col: DG.Column) => ({
        name: col.name,
        diffType: this.activitiesParams[col.name]?.deltaType ?? MmpDiffTypes.delta,
        scaling: this.activitiesParams[col.name]?.scaling ?? SCALING_METHODS.NONE,
      })),
      fragmentCutoff: this.cutoffInput.value,
      runOnFilteredData: this.runOnFilteredDataInput.value,
    });
  }

  loadHistoryString(history: string): void {
    if (!history)
      return;
    try {
      const parsed = JSON.parse(history);
      if (parsed.fragmentCutoff != null)
        this.cutoffInput.value = parsed.fragmentCutoff;
      if (parsed.runOnFilteredData != null)
        this.runOnFilteredDataInput.value = parsed.runOnFilteredData;
      const table = this.tableInput.value;
      if (table != null && Array.isArray(parsed.activities)) {
        // restore only the columns that exist in the current table and are numerical
        const numericalColumns = Array.from(table.columns.numerical);
        const restored: {name: string, diffType: string, scaling: string}[] = [];
        for (const activity of parsed.activities) {
          const col = activity?.name ? table.col(activity.name) : null;
          if (!col || !numericalColumns.includes(col))
            continue;
          const diffType = Object.values(MmpDiffTypes).includes(activity.diffType) ?
            activity.diffType : MmpDiffTypes.delta;
          const allowedScalings = col.stats.min < 0 ? [SCALING_METHODS.NONE] :
            [SCALING_METHODS.NONE, SCALING_METHODS.LG, SCALING_METHODS.MINUS_LG];
          const scaling = allowedScalings.includes(activity.scaling) ? activity.scaling : SCALING_METHODS.NONE;
          restored.push({name: col.name, diffType, scaling});
        }
        // seed the per-activity mappings before setting the selection so the rebuilt choice
        // rows pick them up
        for (const activity of restored)
          this.activitiesParams[activity.name] = {deltaType: activity.diffType, scaling: activity.scaling};
        this.activitiesInput.value = restored.map((activity) => table.col(activity.name));
        // a programmatic value set does not fire onValueChanged — rebuild the rows explicitly
        this.updateActivitiesParamsDiv();
      }
      this.syncCall();
    } catch (e: any) {
      grok.log.error(e);
    }
  }

  inputFor(propertyName: string): DG.InputBase {
    switch (propertyName) {
    case 'table':
      return this.tableInput;
    case 'molecules':
      return this.colInput;
    case 'activities':
      return this.activitiesInput;
    case 'fragmentCutoff':
      return this.cutoffInput;
    case 'runOnFilteredData':
      return this.runOnFilteredDataInput;
    default:
      throw new Error(`Unknown property name: ${propertyName}`);
    }
  }

  get onInputChanged(): Observable<any> {
    return this.inputChangedSubject;
  }
}
