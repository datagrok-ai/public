import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, Subject} from 'rxjs';

import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {scaleActivity} from '../molecular-matched-pairs/mmp-viewer/mmpa-utils';
import {defaultAxis, holdsFragments} from './sar-matrix-columns';
import {MAX_SERIES_LEVELS} from './sar-matrix-run';

const DIRECTIONS = ['Auto (from scaling)', 'Higher is better', 'Lower is better'];
const HISTOGRAM_W = 200;
const HISTOGRAM_H = 110;

/** Dialog of `Chem:SarMatrixAnalysis`: fragment the molecules, or use R-group columns already in the table. */
export class SarMatrixEditor extends DG.FuncCallEditor {
  private readonly tableInput: DG.InputBase<DG.DataFrame | null>;
  private moleculesInput!: DG.InputBase<DG.Column | null>;
  private activityInput!: DG.InputBase<DG.Column | null>;
  private coreInput!: DG.InputBase<DG.Column | null>;
  private rGroupsInput!: DG.InputBase<DG.Column[]>;
  private axisInput!: DG.ChoiceInput<string | null>;
  private seriesInput!: DG.InputBase<DG.Column | null>;
  private readonly moleculesHost = ui.div();
  private readonly activityHost = ui.div();
  private readonly coreHost = ui.div();
  private readonly rGroupsHost = ui.div();
  private readonly axisHost = ui.div();
  private readonly seriesHost = ui.div();
  private readonly layoutHost = ui.divText('', 'chem-sar-editor-layout');
  private readonly rGroupsForm = ui.divV([this.coreHost, this.rGroupsHost, this.axisHost, this.layoutHost]);
  private readonly histogramHost = ui.div([], 'chem-sar-editor-histogram');
  private histogram: DG.Viewer | null = null;

  private readonly scalingInput = ui.input.choice('Scaling', {value: SCALING_METHODS.MINUS_LG,
    items: Object.values(SCALING_METHODS), nullable: false,
    tooltipText: 'Activity transform before the additive model'});
  private readonly directionInput = ui.input.choice('Direction', {value: DIRECTIONS[0], items: DIRECTIONS,
    nullable: false, tooltipText: 'Which end of the scaled activity is more potent'});
  private readonly useRGroupsInput = ui.input.bool('Use existing R-groups', {value: false,
    tooltipText: 'Build the matrices from core and R-group columns already in the table, for example from ' +
      'R-Groups Analysis, instead of fragmenting the molecules',
    onValueChanged: (on: boolean) => this.setMode(on)});
  private readonly cutoffInput = ui.input.float('Fragment cutoff', {value: 0.4, min: 0.1, max: 1,
    nullable: false, tooltipText: 'Largest substituent kept when fragmenting, as a fraction of the core'});
  private readonly levelsInput = ui.input.int('Series levels', {value: 3, min: 1, max: MAX_SERIES_LEVELS,
    nullable: false, tooltipText: 'Nested matrix tiers (L1, L2, …); each level folds matrices one cut broader'});
  private readonly mcsInput = ui.input.bool('Group leftovers by MCS', {value: false,
    tooltipText: 'Also cover the compounds no shared core could group, by searching those for a common core'});
  private readonly predictInput = ui.input.bool('Predict analogs', {value: true,
    tooltipText: 'Fill unmade combinations with Free-Wilson predictions'});
  private readonly fragmentationForm: HTMLElement;
  private readonly settingsIcon: HTMLElement;
  /** Inputs that map one-to-one onto a function parameter. */
  private readonly simple: [string, DG.InputBase<any>][] = [['scaling', this.scalingInput],
    ['activityDirection', this.directionInput], ['fragmentCutoff', this.cutoffInput],
    ['fragmentationLevels', this.levelsInput], ['predictVirtual', this.predictInput],
    ['useMcsAnchors', this.mcsInput]];

  private readonly inputChanged = new Subject<any>();

  constructor(private readonly funcCall: DG.FuncCall) {
    const root = ui.div([]);
    super(root);
    this.fragmentationForm = ui.form([this.cutoffInput, this.levelsInput, this.mcsInput]);
    this.fragmentationForm.style.display = 'none';
    this.settingsIcon = ui.icons.settings(() => {
      const shown = this.fragmentationForm.style.display !== 'none';
      this.fragmentationForm.style.display = shown ? 'none' : 'flex';
      this.fragmentationForm.classList.remove('ui-form-condensed');
    }, 'Fragmentation settings');
    this.useRGroupsInput.root.appendChild(this.settingsIcon);
    this.tableInput = ui.input.table('Table', {
      value: funcCall.inputs['table'] ?? grok.shell.tv?.dataFrame,
      items: grok.shell.tables,
      onValueChanged: () => this.onTableChanged(),
    });

    for (const [, input] of this.simple)
      input.onChanged.subscribe(() => this.syncCall());
    this.scalingInput.onChanged.subscribe(() => this.drawHistogram());

    this.rGroupsForm.style.display = 'none';
    this.onTableChanged();
    root.append(this.getEditor());
  }

  /** Column inputs are bound to one table. */
  private onTableChanged(): void {
    const table = this.tableInput.value;
    if (table === null)
      return;
    this.moleculesInput = ui.input.column('Molecules', {table, nullable: false,
      value: table.columns.toList().find((c) => c.semType === DG.SEMTYPE.MOLECULE && !holdsFragments(c)),
      filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
      onValueChanged: () => this.syncCall()});
    this.activityInput = ui.input.column('Activity', {table, nullable: false,
      value: table.columns.toList().find((c) => this.isActivity(c)),
      filter: (c: DG.Column) => this.isActivity(c),
      onValueChanged: () => {
        this.onActivityChanged();
        this.syncCall();
      }});
    this.coreInput = ui.input.column('Core', {table, filter: holdsFragments,
      tooltipText: 'Column with the core, its attachment points marked [*:1], [*:2], …',
      onValueChanged: () => this.onRGroupsChanged()});
    this.rGroupsInput = ui.input.columns('R-groups', {table, value: [], filter: holdsFragments,
      tooltipText: 'Columns with the substituent at each attachment point of the core',
      onValueChanged: () => this.onRGroupsChanged()});
    this.axisInput = ui.input.choice('Matrix columns', {value: null, items: [],
      tooltipText: 'The R-group whose substituents become the matrix columns. ' +
        'The core and the other R-groups make up the rows',
      onValueChanged: () => {
        this.describeLayout();
        this.syncCall();
      }});
    this.seriesInput = ui.input.column('Series', {table,
      tooltipText: 'Optional. Compounds sharing a value make one matrix, named with that value',
      onValueChanged: () => {
        this.describeLayout();
        this.syncCall();
      }});
    // A column input starts on the table's first matching column; these start empty.
    for (const input of [this.coreInput, this.seriesInput]) {
      input.nullable = true;
      input.value = null;
    }

    for (const [host, input] of [[this.moleculesHost, this.moleculesInput], [this.activityHost, this.activityInput],
      [this.coreHost, this.coreInput], [this.rGroupsHost, this.rGroupsInput],
      [this.axisHost, this.axisInput], [this.seriesHost, this.seriesInput]] as
      [HTMLElement, DG.InputBase][]) {
      ui.empty(host);
      host.appendChild(input.root);
    }
    this.onActivityChanged();
    this.onRGroupsChanged();
  }

  /** A DateTime column is numerical, but potency arithmetic on a timestamp is nonsense. */
  private isActivity(column: DG.Column): boolean {
    return column.isNumerical && column.type !== DG.COLUMN_TYPE.DATE_TIME;
  }

  /** A log scale is offered only for positive values. */
  private onActivityChanged(): void {
    const column = this.activityInput?.value;
    const scalable = column !== null && column !== undefined && column.stats.min > 0;
    this.scalingInput.enabled = scalable;
    if (!scalable)
      this.scalingInput.value = SCALING_METHODS.NONE;
    this.drawHistogram();
  }

  private drawHistogram(): void {
    const column = this.activityInput?.value;
    this.histogram?.detach();
    this.histogram = null;
    ui.empty(this.histogramHost);
    if (!column)
      return;
    const scaled = scaleActivity(column as DG.Column<number>, this.scalingInput.value ?? SCALING_METHODS.NONE);
    this.histogram = DG.DataFrame.fromColumns([scaled]).plot.histogram({
      filteringEnabled: false, legendVisibility: 'Never', showXAxis: true,
      showColumnSelector: false, showRangeSlider: false, showBinSelector: false,
    });
    this.histogram.root.style.width = `${HISTOGRAM_W}px`;
    this.histogram.root.style.height = `${HISTOGRAM_H}px`;
    this.histogramHost.appendChild(this.histogram.root);
  }

  private rGroupColumns(): DG.Column[] {
    return this.rGroupsInput?.value ?? [];
  }

  private usesRGroups(): boolean {
    return this.useRGroupsInput.value;
  }

  /** The picked columns are kept while hidden, so switching back restores them. */
  private setMode(rGroups: boolean): void {
    this.rGroupsForm.style.display = rGroups ? 'flex' : 'none';
    this.settingsIcon.style.display = rGroups ? 'none' : '';
    if (rGroups)
      this.fragmentationForm.style.display = 'none';
    this.onRGroupsChanged();
  }

  /** Matrix columns offers the picked R-groups; a pick survives as long as its column stays picked. */
  private onRGroupsChanged(): void {
    if (this.axisInput === undefined)
      return;
    const columns = this.rGroupColumns();
    const names = columns.map((c) => c.name);
    const chosen = this.axisInput.value;
    this.axisInput.items = names;
    if (names.length === 0)
      this.axisInput.value = null;
    else if (chosen === null || !names.includes(chosen))
      this.axisInput.value = defaultAxis(columns)?.name ?? names[names.length - 1];
    this.describeLayout();
    this.syncCall();
  }

  private describeLayout(): void {
    const core = this.coreInput?.value ?? null;
    const axis = this.axisInput?.value ?? null;
    if (core === null || axis === null) {
      this.layoutHost.textContent = core === null ? 'Pick the core column.' : 'Pick the R-groups.';
      return;
    }
    const rows = [core.name, ...this.rGroupColumns().map((c) => c.name).filter((name) => name !== axis)];
    const series = this.seriesInput?.value?.name;
    this.layoutHost.textContent = `${series ? `One matrix per ${series}` : 'One matrix per core'}  ·  ` +
      `Rows: ${rows.join(' + ')}  ·  Columns: ${axis}`;
  }

  private syncCall(): void {
    const rGroups = this.usesRGroups();
    this.funcCall.inputs['table'] = this.tableInput.value;
    this.funcCall.inputs['molecules'] = this.moleculesInput?.value;
    this.funcCall.inputs['activity'] = this.activityInput?.value;
    for (const [key, input] of this.simple)
      this.funcCall.inputs[key] = input.value;
    this.funcCall.inputs['seriesColumn'] = this.seriesInput?.value?.name ?? '';
    this.funcCall.inputs['coreColumn'] = rGroups ? this.coreInput?.value ?? null : null;
    this.funcCall.inputs['rGroupColumns'] = rGroups ? this.rGroupColumns() : [];
    this.funcCall.inputs['matrixColumns'] = rGroups ? this.axisInput?.value ?? '' : '';
    this.inputChanged.next(null);
  }

  getEditor(): HTMLElement {
    return ui.divV([
      ui.divH([
        ui.divV([this.tableInput.root, this.moleculesHost, this.activityHost, this.scalingInput.root,
          this.directionInput.root], {style: {flex: '1 1 auto'}}),
        this.histogramHost,
      ]),
      this.useRGroupsInput.root,
      this.fragmentationForm,
      this.rGroupsForm,
      this.seriesHost,
      this.predictInput.root,
    ], {style: {minWidth: '440px'}});
  }

  get isValid(): boolean {
    const core = this.coreInput?.value ?? null;
    const axis = this.axisInput?.value ?? null;
    if (this.usesRGroups() && (core === null || axis === null || core.name === axis))
      return false;
    // A pure read: the platform re-evaluates it on every onInputChanged emission.
    return this.tableInput.value !== null && (this.moleculesInput?.value ?? null) !== null &&
      (this.activityInput?.value ?? null) !== null && this.cutoffInput.value !== null &&
      this.levelsInput.value !== null;
  }

  getHistoryString(): string {
    return JSON.stringify({
      ...Object.fromEntries(this.simple.map(([key, input]) => [key, input.value])),
      series: this.seriesInput?.value?.name ?? null,
      core: this.coreInput?.value?.name ?? null,
      rgroups: this.rGroupColumns().map((c) => c.name),
      axis: this.axisInput?.value ?? null,
    });
  }

  loadHistoryString(history: string): void {
    if (!history)
      return;
    try {
      const parsed = JSON.parse(history);
      for (const [key, input] of this.simple) {
        if (parsed[key] != null)
          input.value = parsed[key];
      }
      const table = this.tableInput.value;
      if (table !== null) {
        this.seriesInput.value = parsed.series == null ? null : table.col(parsed.series);
        const core = parsed.core == null ? null : table.col(parsed.core);
        const rgroups = (parsed.rgroups ?? [])
          .map((name: string) => table.col(name))
          .filter((c: DG.Column | null) => c !== null);
        this.useRGroupsInput.value = core !== null || rgroups.length > 0;
        this.coreInput.value = core;
        this.rGroupsInput.value = rgroups;
        // The choice accepts the axis only once it offers that name.
        this.onRGroupsChanged();
        if (parsed.axis != null && rgroups.some((c: DG.Column) => c.name === parsed.axis))
          this.axisInput.value = parsed.axis;
      }
      this.syncCall();
    } catch (e: any) {
      grok.log.error(e);
    }
  }

  inputFor(propertyName: string): DG.InputBase {
    const simple = this.simple.find(([key]) => key === propertyName);
    if (simple !== undefined)
      return simple[1];
    switch (propertyName) {
    case 'table':
      return this.tableInput;
    case 'molecules':
      return this.moleculesInput;
    case 'activity':
      return this.activityInput;
    case 'coreColumn':
      return this.coreInput;
    case 'rGroupColumns':
      return this.rGroupsInput;
    case 'matrixColumns':
      return this.axisInput;
    case 'seriesColumn':
      return this.seriesInput;
    default:
      throw new Error(`Unknown property name: ${propertyName}`);
    }
  }

  get onInputChanged(): Observable<any> {
    return this.inputChanged;
  }
}
