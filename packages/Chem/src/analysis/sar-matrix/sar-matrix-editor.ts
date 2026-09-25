import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, Subject} from 'rxjs';

import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {scaleActivity} from '../molecular-matched-pairs/mmp-viewer/mmpa-utils';
import {defaultAxis} from './sar-matrix-columns';
import {MAX_SERIES_LEVELS} from './sar-matrix-run';

const DIRECTIONS = ['Auto (from scaling)', 'Higher is better', 'Lower is better'];

/**
 * FuncCall editor for `Chem:SarMatrixAnalysis`.
 *
 * Two ways to get the series, behind one checkbox: fragment the structures, or read a decomposition
 * the table already holds. The decomposition is named the way it is written — the scaffold, then its
 * R-groups, then which R-group the matrix enumerates across — and naming one hides the fragmentation
 * settings, which describe a split no longer being made.
 */
export class SarMatrixEditor extends DG.FuncCallEditor {
  private readonly tableInput: DG.InputBase<DG.DataFrame | null>;
  private moleculesInput!: DG.InputBase<DG.Column | null>;
  private activityInput!: DG.InputBase<DG.Column | null>;
  private coreInput!: DG.InputBase<DG.Column | null>;
  private fragmentsInput!: DG.InputBase<DG.Column[]>;
  private axisInput!: DG.ChoiceInput<string | null>;
  private seriesInput!: DG.InputBase<DG.Column | null>;
  private readonly moleculesHost = ui.div();
  private readonly activityHost = ui.div();
  private readonly coreHost = ui.div();
  private readonly fragmentsHost = ui.div();
  private readonly axisHost = ui.div();
  private readonly seriesHost = ui.div();
  private readonly layoutHost = ui.divText('', 'chem-sar-editor-layout');
  /** The three pickers, revealed together: naming a core without its R-groups describes nothing. */
  private readonly decompositionForm = ui.divV([this.coreHost, this.fragmentsHost, this.axisHost,
    this.layoutHost]);
  private readonly histogramHost = ui.div([], 'chem-sar-editor-histogram');

  private readonly scalingInput = ui.input.choice('Scaling', {value: SCALING_METHODS.MINUS_LG,
    items: Object.values(SCALING_METHODS), nullable: false,
    tooltipText: 'Activity transform before the additive model'});
  private readonly directionInput = ui.input.choice('Direction', {value: DIRECTIONS[0], items: DIRECTIONS,
    nullable: false, tooltipText: 'Which end of the scaled activity is more potent'});
  private readonly fragmentStructuresInput = ui.input.bool('Fragment structures', {value: true,
    tooltipText: 'Cut the molecules into cores and substituents. Turn off to name columns that ' +
      'already hold a decomposition',
    onValueChanged: (on: boolean) => this.setMode(!on)});
  private readonly cutoffInput = ui.input.float('Fragment cutoff', {value: 0.4, min: 0.1, max: 1,
    nullable: false, tooltipText: 'Largest substituent kept when fragmenting, as a fraction of the core'});
  private readonly levelsInput = ui.input.int('Series levels', {value: 3, min: 1, max: MAX_SERIES_LEVELS,
    nullable: false, tooltipText: 'Nested matrix tiers (L1, L2, …); each level folds matrices one cut broader'});
  private readonly mcsInput = ui.input.bool('Group leftovers by MCS', {value: false,
    tooltipText: 'Also cover the compounds no shared core could group, by searching those for a common core'});
  private readonly predictInput = ui.input.bool('Predict analogs', {value: true,
    tooltipText: 'Fill unmade row × column cells with Free-Wilson predictions'});
  private readonly fragmentationForm: HTMLElement;
  private readonly settingsIcon: HTMLElement;
  /** Inputs that map one-to-one onto a function parameter, so they sync, save and resolve alike. */
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
    this.fragmentStructuresInput.root.appendChild(this.settingsIcon);
    this.tableInput = ui.input.table('Table', {
      value: funcCall.inputs['table'] ?? grok.shell.tv?.dataFrame,
      items: grok.shell.tables,
      onValueChanged: () => this.onTableChanged(),
    });

    for (const [, input] of this.simple)
      input.onChanged.subscribe(() => this.syncCall());
    this.scalingInput.onChanged.subscribe(() => this.drawHistogram());

    this.decompositionForm.style.display = 'none';
    this.onTableChanged();
    root.append(this.getEditor());
  }

  /** Column inputs are bound to one table, so a different table means new ones. */
  private onTableChanged(): void {
    const table = this.tableInput.value;
    if (table === null)
      return;
    this.moleculesInput = ui.input.column('Molecules', {table, nullable: false,
      value: table.columns.toList().find((c) => c.semType === DG.SEMTYPE.MOLECULE),
      filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
      onValueChanged: () => this.syncCall()});
    this.activityInput = ui.input.column('Activity', {table, nullable: false,
      value: table.columns.toList().find((c) => this.isActivity(c)),
      filter: (c: DG.Column) => this.isActivity(c),
      onValueChanged: () => {
        this.onActivityChanged();
        this.syncCall();
      }});
    this.coreInput = ui.input.column('Core column', {table,
      tooltipText: 'The scaffold every row is drawn from',
      onValueChanged: () => this.onFragmentsChanged()});
    this.coreInput.nullable = true;
    this.coreInput.value = null;
    this.fragmentsInput = ui.input.columns('R-group columns', {table, value: [],
      tooltipText: 'The substituent columns that hang off the core — R1, R2, … or whatever they are called',
      onValueChanged: () => this.onFragmentsChanged()});
    this.axisInput = ui.input.choice('Columns axis', {value: null, items: [],
      tooltipText: 'The R-group the matrix enumerates across. The rest fold into the row',
      onValueChanged: () => {
        this.describeLayout();
        this.syncCall();
      }});
    this.seriesInput = ui.input.column('Series column', {table,
      tooltipText: 'Optional. Compounds sharing a value make one matrix, named with that value',
      onValueChanged: () => {
        this.describeLayout();
        this.syncCall();
      }});
    // A column input resolves to the first column of the table on construction, which would group
    // every compound by something nobody asked for. Emptied after the fact, as nullable alone does not.
    this.seriesInput.nullable = true;
    this.seriesInput.value = null;

    for (const [host, input] of [[this.moleculesHost, this.moleculesInput], [this.activityHost, this.activityInput],
      [this.coreHost, this.coreInput], [this.fragmentsHost, this.fragmentsInput],
      [this.axisHost, this.axisInput], [this.seriesHost, this.seriesInput]] as
      [HTMLElement, DG.InputBase][]) {
      ui.empty(host);
      host.appendChild(input.root);
    }
    this.onActivityChanged();
    this.onFragmentsChanged();
  }

  /** A DateTime column reports isNumerical — dates are numeric internally, which is what lets them
   *  serve as a plot axis — and potency arithmetic on a timestamp is silent nonsense. */
  private isActivity(column: DG.Column): boolean {
    return column.isNumerical && column.type !== DG.COLUMN_TYPE.DATE_TIME;
  }

  /** A log scale turns a zero or a negative into ±Infinity, so it is offered only where it applies. */
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
    ui.empty(this.histogramHost);
    if (!column)
      return;
    const scaled = scaleActivity(column as DG.Column<number>, this.scalingInput.value ?? SCALING_METHODS.NONE);
    const histogram = DG.DataFrame.fromColumns([scaled]).plot.histogram({
      filteringEnabled: false, legendVisibility: 'Never', showXAxis: true,
      showColumnSelector: false, showRangeSlider: false, showBinSelector: false,
    });
    // A viewer lays out at its own default height, not its host's, and spills over the inputs below
    // it: the host is a fixed strip beside the form, not a pane the viewer can size itself in.
    histogram.root.style.height = '100%';
    this.histogramHost.appendChild(histogram.root);
  }

  private fragmentColumns(): DG.Column[] {
    return this.fragmentsInput?.value ?? [];
  }

  /** Whether the table's own decomposition is what builds the series. */
  private usesColumns(): boolean {
    return !this.fragmentStructuresInput.value;
  }

  /** What was picked is kept while hidden — `syncCall` ignores it — so turning fragmentation back off
   *  restores the pairing instead of making the user name it again. */
  private setMode(columns: boolean): void {
    this.decompositionForm.style.display = columns ? 'flex' : 'none';
    // Neither the fragmentation settings nor the gear that reopens them describe what will be run.
    this.settingsIcon.style.display = columns ? 'none' : '';
    if (columns)
      this.fragmentationForm.style.display = 'none';
    this.onFragmentsChanged();
  }

/** The axis choice offers exactly the R-groups currently picked, defaulting per {@link defaultAxis} to
 *  the layout an R-group decomposition asks for. A pick the user made by hand survives as long as its
 *  column does. */
  private onFragmentsChanged(): void {
    if (this.axisInput === undefined)
      return;
    const columns = this.fragmentColumns();
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

  /** States the layout in words rather than a count: what goes down the side and what goes across is
   *  the question the three pickers raise, and it is answerable without touching the data. */
  private describeLayout(): void {
    const core = this.coreInput?.value ?? null;
    const columns = this.fragmentColumns();
    const axis = this.axisInput?.value ?? null;
    if (core === null || axis === null) {
      this.layoutHost.textContent = core === null ?
        'Name the core column — the scaffold every row is drawn from.' :
        'Name the R-group the matrix enumerates across.';
      return;
    }
    const rows = columns.filter((c) => c.name !== axis).map((c) => c.name);
    const seriesBy = this.seriesInput?.value?.name;
    this.layoutHost.textContent = seriesBy ?
      `One series per ${seriesBy}  ·  Rows: ${[core.name, ...rows].join(' + ')}  ·  Columns: ${axis}` :
      rows.length > 0 ?
        `One series per ${core.name}  ·  Rows: ${rows.join(' + ')}  ·  Columns: ${axis}` :
        `Rows: ${core.name}  ·  Columns: ${axis}`;
  }

  private syncCall(): void {
    const columns = this.usesColumns();
    this.funcCall.inputs['table'] = this.tableInput.value;
    this.funcCall.inputs['molecules'] = this.moleculesInput?.value;
    this.funcCall.inputs['activity'] = this.activityInput?.value;
    for (const [key, input] of this.simple)
      this.funcCall.inputs[key] = input.value;
    this.funcCall.inputs['seriesColumn'] = this.seriesInput?.value?.name ?? '';
    this.funcCall.inputs['coreColumn'] = columns ? this.coreInput?.value ?? null : null;
    this.funcCall.inputs['fragmentColumns'] = columns ? this.fragmentColumns() : [];
    this.funcCall.inputs['columnAxis'] = columns ? this.axisInput?.value ?? '' : '';
    this.inputChanged.next(null);
  }

  getEditor(): HTMLElement {
    return ui.divV([
      ui.divH([
        ui.divV([this.tableInput.root, this.moleculesHost, this.activityHost, this.scalingInput.root,
          this.directionInput.root], {style: {flex: '1 1 auto'}}),
        this.histogramHost,
      ]),
      this.fragmentStructuresInput.root,
      this.fragmentationForm,
      this.decompositionForm,
      this.seriesHost,
      this.predictInput.root,
    ], {style: {minWidth: '440px'}});
  }

  get isValid(): boolean {
    // Half a decomposition builds nothing: the core says what the rows share, the axis what they
    // vary, and one column cannot be both.
    const core = this.coreInput?.value ?? null;
    const axis = this.axisInput?.value ?? null;
    if (this.usesColumns() && (core === null || axis === null || core.name === axis))
      return false;
    // Pure read: the platform re-evaluates this on every onInputChanged emission, so emitting from
    // here would spin. The handlers keep the inputs current.
    return this.tableInput.value !== null && (this.moleculesInput?.value ?? null) !== null &&
      (this.activityInput?.value ?? null) !== null && this.cutoffInput.value !== null &&
      this.levelsInput.value !== null;
  }

  getHistoryString(): string {
    return JSON.stringify({
      ...Object.fromEntries(this.simple.map(([key, input]) => [key, input.value])),
      core: this.coreInput?.value?.name ?? null,
      rgroups: this.fragmentColumns().map((c) => c.name),
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
        const core = parsed.core == null ? null : table.col(parsed.core);
        const rgroups = (parsed.rgroups ?? [])
          .map((name: string) => table.col(name))
          .filter((c: DG.Column | null) => c !== null);
        // Set both ways: a run that fragmented must reopen as one, whatever mode the dialog is in.
        this.fragmentStructuresInput.value = core === null && rgroups.length === 0;
        this.coreInput.value = core;
        this.fragmentsInput.value = rgroups;
        // Before the axis, which the choice only accepts once it offers that name.
        this.onFragmentsChanged();
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
    case 'fragmentColumns':
      return this.fragmentsInput;
    case 'columnAxis':
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
