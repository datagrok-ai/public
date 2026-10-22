import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {signal, div, divV, span, text, plural, badge, notify, Filters, Wizard, DataTable, ChoiceInput, Form} from
  '@datagrok-libraries/u2';
import type {BadgeVariant, DataTableColumn, IProperty} from '@datagrok-libraries/u2';
import {fromDartInput, moleculeRenderer} from '@datagrok-libraries/u2/src/dg/index.js';
import {DomainErrors, DomainPick, DomainTable, ImportMapping} from '@datagrok-libraries/ems/src/dg';
import {MolTrackEms} from './app';
import {Registration} from './registration';
import type {RegistrationMapping, RegistrationPreview, RegistrationStatus} from './registration';
import {MolTrackTables} from './tables';

export type RegisterTarget = 'Compounds' | 'Batches';

const STRUCTURE = 'Structure';
const SKIP_ROW = 'Skip the row';
const REJECT_ALL = 'Reject all';
const VERDICTS: Record<RegistrationPreview['verdict'], {label: string, variant: BadgeVariant}> = {
  new: {label: 'New', variant: 'accent'},
  existing: {label: 'Existing', variant: 'default'},
  error: {label: 'Error', variant: 'error'},
};
const OUTCOMES: Record<RegistrationStatus, {label: string, variant: BadgeVariant}> = {
  success: {label: 'Registered', variant: 'success'},
  failed: {label: 'Failed', variant: 'error'},
  not_processed: {label: 'Not registered', variant: 'default'},
};
const THUMBNAIL = moleculeRenderer({width: 80, height: 40});

/** Bulk registration as the legacy app offered it: a table or a file, its columns mapped onto the compound and
 * the batch, a preview of what each row becomes, and the per-row report. */
export class RegisterDialog {
  private readonly _frame = signal<DG.DataFrame | null>(null);
  private readonly _columns = new ImportMapping();
  /** What each target reads in its choice: the gate names a target claimed twice by it. */
  private _labels = new Map<string, string>();
  private _tables: {compounds: DomainTable, batches?: DomainTable} | undefined;
  private readonly _preview = div([], 'moltrack-ems-register-preview');
  private readonly _report = div([], 'moltrack-ems-register-report');
  private readonly _wizard: Wizard;
  private readonly _target: ChoiceInput;
  private readonly _onError: ChoiceInput;
  private _previewTable: DataTable<number> | undefined;
  private _result: DG.DataFrame | null = null;
  private _remaps = 0;
  private _previews = 0;
  private _settled = false;

  private constructor(source: DG.DataFrame | null, target: RegisterTarget,
    private readonly _resolve: (report: DG.DataFrame | null) => void) {
    this._frame.value = source;
    // the wizard builds its first step at once, and that step holds these
    this._target = new ChoiceInput({label: 'Register', value: target, nullable: false, items: ['Compounds', 'Batches'],
      onChanged: () => void this._remap()});
    // a run that registered some rows cannot be run again as it is: the rows it wrote would be registered twice
    this._onError = new ChoiceInput({label: 'On error', value: REJECT_ALL, nullable: false,
      items: [REJECT_ALL, SKIP_ROW],
      tooltipText: `${REJECT_ALL} registers nothing while a row has a problem; ${SKIP_ROW} registers the others`});
    this._wizard = new Wizard({
      steps: [
        {id: 'source', title: 'Source', content: () => this._sourceStep(),
          canProceed: () => this._frame.value === null ? 'Choose an open table or a file' : null},
        {id: 'mapping', title: 'Mapping', content: () => this._mappingStep(), canProceed: () => this._problem()},
        {id: 'preview', title: 'Preview', content: this._preview, canProceed: () => this._problem(),
          nextText: 'Register', onActivate: () => void this._showPreview()},
        {id: 'report', title: 'Report', content: this._report, done: true, onActivate: () => void this._register()},
      ],
      onFinish: () => this._settle(this._result),
      onCancel: () => this._settle(null),
    });
    this._wizard.root.classList.add('moltrack-ems-register');
  }

  /** Resolves to the registration report, or null when cancelled. */
  static open(source: DG.DataFrame | null = grok.shell.t ?? grok.shell.tables[0] ?? null,
    target: RegisterTarget = 'Batches'): Promise<DG.DataFrame | null> {
    return new Promise((resolve) => {
      const dialog = new RegisterDialog(source, target, resolve);
      dialog._remap().then(() => dialog._wizard.openInDialog('Register in bulk', {width: 860, height: 600}), (e) => {
        notify.error(DomainErrors.message(e));
        dialog._settle(null);
      });
    });
  }

  private _settle(report: DG.DataFrame | null): void {
    if (this._settled)
      return;
    this._settled = true;
    this._columns.dispose();
    this._target.dispose();
    this._onError.dispose();
    this._wizard.dispose();
    this._resolve(report);
  }

  private _sourceStep(): HTMLElement {
    const picker = fromDartInput<DG.DataFrame | null>(ui.input.table('Data', {value: this._frame.peek() ?? undefined}),
      'source');
    picker.root.title = 'Choose an open table or a file';
    const summary = span('', 'moltrack-ems-register-summary');
    const form = new Form().addAll([picker, this._target, this._onError]);
    form.effect(() => {
      const frame = picker.value.value;
      summary.textContent = frame === null ? '' :
        `${frame.name}: ${plural(frame.rowCount, 'row', 'rows')}, ${plural(frame.columns.length, 'column', 'columns')}`;
      if (frame !== this._frame.peek()) {
        this._frame.value = frame;
        void this._remap();
      }
    });
    return divV([form.root, summary], 'moltrack-ems-register-source');
  }

  private _mappingStep(): HTMLElement {
    const table = this._columns.table('Column in the file', 'Registers as');
    table.root.classList.add('moltrack-ems-register-mapping');
    return table.root;
  }

  /** Again on every visit: going back may have changed the file, the target, the mapping or the On error. */
  private async _showPreview(): Promise<void> {
    const at = ++this._previews;
    const frame = this._frame.peek()!;
    const mapping = this._mapping();
    this._preview.replaceChildren(span('Checking the structures…', 'moltrack-ems-register-summary'));
    try {
      const rows = await Registration.preview(frame, mapping);
      if (at !== this._previews || this._settled)
        return;
      const structure = frame.col(mapping.structure)!;
      this._previewTable?.dispose();
      const table = this._previewTable = this._wizard.runInScope(() => new DataTable<number>({
        rowHeight: 44,
        keyOf: (i) => String(i),
        columns: [
          RegisterDialog._rowNumber(),
          {name: 'structure', header: 'Structure', width: '96px', render: (i) => rows[i].verdict === 'error' ?
            text(structure.get(i)) : THUMBNAIL.listItem!(structure.get(i))},
          ...RegisterDialog._fileColumns(frame, this._identifying(mapping)),
          {name: 'verdict', header: 'Becomes', width: '84px', render: (i) =>
            badge(VERDICTS[rows[i].verdict].label, {variant: VERDICTS[rows[i].verdict].variant})},
          {name: 'message', header: '', width: '2fr', render: (i) => rows[i].message},
        ],
      }));
      table.setItems(Array.from({length: rows.length}, (_, i) => i));
      table.root.classList.add('moltrack-ems-register-rows');
      this._preview.replaceChildren(...this._summary(rows, mapping.batch !== undefined), table.root);
    } catch (e) {
      if (at === this._previews)
        this._preview.replaceChildren(span(DomainErrors.message(e), 'moltrack-ems-register-problem'));
    }
  }

  /** What the rows become, then what On error does with the rows that have a problem. */
  private _summary(rows: RegistrationPreview[], batches: boolean): HTMLElement[] {
    const count = (test: (row: RegistrationPreview) => boolean) => rows.filter(test).length;
    const compounds = count((row) => row.verdict === 'new' && row.sameAs === undefined);
    const problems = count((row) => row.verdict === 'error');
    const registered = count((row) => row.verdict === 'existing');
    const parts = [plural(compounds, 'new compound', 'new compounds')];
    if (batches) {
      parts.push(plural(rows.length - problems, 'batch', 'batches') + (registered === 0 ? '' :
        ` (${registered} of ${registered === 1 ? 'a registered compound' : 'registered compounds'})`));
    }
    if (problems > 0)
      parts.push(`${problems} with a problem`);
    const lines = [span(`${plural(rows.length, 'row', 'rows')}: ${parts.join(', ')}`, 'moltrack-ems-register-summary')];
    const rejectAll = this._onError.value.peek() === REJECT_ALL;
    if (problems > 0) {
      lines.push(rejectAll ? span(`Nothing is registered while a row has a problem: fix the file, or choose ` +
        `"${SKIP_ROW}" on the first step.`, 'moltrack-ems-register-problem') :
        span('The rows with a problem are skipped.', 'moltrack-ems-register-summary'));
    }
    const writes = compounds + (batches ? rows.length : 0);
    if (rejectAll && problems === 0 && writes > Registration.MAX_OPS) {
      lines.push(span(`${REJECT_ALL} writes at most ${Registration.MAX_OPS} rows at once and this file needs ` +
        `${writes}: split the file, or choose "${SKIP_ROW}" on the first step.`, 'moltrack-ems-register-problem'));
    }
    return lines;
  }

  private async _register(): Promise<void> {
    const frame = this._frame.peek()!;
    const mapping = this._mapping();
    this._report.replaceChildren(span('Registering…', 'moltrack-ems-register-summary'));
    try {
      const report = this._result = await Registration.register(frame, mapping,
        {rejectAll: this._onError.value.peek() === REJECT_ALL});
      if (this._settled)
        return;
      report.name = `${frame.name} registered`;
      const status = report.col(Registration.STATUS)!;
      const ok = status.toList().filter((s) => s === 'success').length;
      const table = this._wizard.runInScope(() => new DataTable<number>({
        rowHeight: 24,
        keyOf: (i) => String(i),
        columns: [
          RegisterDialog._rowNumber(),
          ...RegisterDialog._fileColumns(report, this._identifying(mapping)),
          {name: Registration.STATUS, header: 'Result', width: '120px', render: (i) =>
            badge(OUTCOMES[status.get(i) as RegistrationStatus].label,
              {variant: OUTCOMES[status.get(i) as RegistrationStatus].variant})},
          ...[Registration.COMPOUND_ID, Registration.BATCH_ID].filter((name) => report.col(name) !== null)
            .map((name) => ({name, header: name, width: '110px', render: (i: number) => text(report.get(name, i))})),
          {name: Registration.MESSAGE, header: '', width: '2fr',
            render: (i) => text(report.get(Registration.MESSAGE, i))},
        ],
      }));
      table.setItems(Array.from({length: report.rowCount}, (_, i) => i));
      table.root.classList.add('moltrack-ems-register-rows');
      const registered = RegisterDialog.registered(report);
      this._report.replaceChildren(span(`Registered ${ok} of ${plural(report.rowCount, 'row', 'rows')}`,
        'moltrack-ems-register-summary'), table.root,
      ...(registered === null ? [] : [ui.link(`Show the ${plural(registered.count, ...registered.noun)} registered`,
        () => this._show(registered.segment, registered.query))]),
      ui.link('Open the report as a table', () => grok.shell.addTableView(report)));
    } catch (e) {
      this._report.replaceChildren(span(DomainErrors.message(e), 'moltrack-ems-register-problem'));
    }
  }

  /** The list that holds what a run registered — its batches, or its compounds — as the view and the filter that
   * show it; null when nothing was registered. */
  static registered(report: DG.DataFrame):
    {segment: string, query: string, count: number, noun: [string, string]} | null {
    const batches = report.col(Registration.BATCH_ID) !== null;
    const ids = report.col(batches ? Registration.BATCH_ID : Registration.COMPOUND_ID)!;
    const status = report.col(Registration.STATUS)!;
    const numbers = [...new Set(Array.from({length: report.rowCount}, (_, i) => i)
      .filter((i) => status.get(i) === 'success' && ids.get(i) !== '')
      .map((i) => Number(DomainPick.digits(ids.get(i)))))];
    if (numbers.length === 0)
      return null;
    return {segment: batches ? 'Batches' : 'Compounds', count: numbers.length,
      noun: batches ? ['batch', 'batches'] : ['compound', 'compounds'],
      query: Filters.format(Filters.cond(batches ? 'batch_regno' : 'molregno', 'in', numbers))};
  }

  /** Closes the report and filters the list of `segment` by `query`. */
  private _show(segment: string, query: string): void {
    this._wizard.next();
    void MolTrackEms.show(`${segment}?q=${encodeURIComponent(query)}`);
  }

  /** The file's columns a person tells its rows apart by: the first one mapped onto a column the compound is
   * searched by, and the first one onto a column the batch is. */
  private _identifying(mapping: RegistrationMapping): string[] {
    const first = (fields: Record<string, string> | undefined, table: DomainTable | undefined) =>
      Object.entries(fields ?? {}).find(([, target]) => table?.info.searchableColumns.includes(target))?.[0];
    return [first(mapping.compound, this._tables?.compounds), first(mapping.batch, this._tables?.batches)]
      .filter((name): name is string => name !== undefined);
  }

  private static _rowNumber(): DataTableColumn<number> {
    return {name: 'row', header: '#', width: '40px', render: (i) => String(i + 1)};
  }

  private static _fileColumns(frame: DG.DataFrame, names: string[]): DataTableColumn<number>[] {
    return names.map((name) => ({name, header: name, render: (i: number) =>
      text(DomainTable.cellValue(frame.col(name)!, i))}));
  }

  /** The choices of every source column, auto-matched by name and caption. */
  private async _remap(): Promise<void> {
    const at = ++this._remaps;
    const frame = this._frame.peek();
    const batches = this._target.value.peek() === 'Batches';
    const compounds = await MolTrackTables.compounds();
    const batchTable = batches ? await MolTrackTables.batches() : undefined;
    if (at !== this._remaps || this._settled)
      return;
    const option = (prefix: string, p: IProperty) => ({value: `${prefix}.${p.name}`,
      label: `${prefix === 'compound' ? 'Compound' : 'Batch'} · ${p.friendlyName ?? p.name}`});
    const targets = [{value: STRUCTURE, label: STRUCTURE},
      ...Registration.targets(compounds).map((p) => option('compound', p)),
      ...(batchTable === undefined ? [] : Registration.targets(batchTable).map((p) => option('batch', p)))];
    const auto = frame === null ? null : Registration.autoMap(frame, compounds, batchTable);
    this._tables = {compounds, batches: batchTable};
    this._labels = new Map(targets.map((t) => [t.value, t.label]));
    this._columns.remap(frame?.columns.names().filter((n) => !Registration.OUTCOME.includes(n)) ?? [], targets,
      (name: string) => auto!.structure === name ? STRUCTURE : auto!.batch?.[name] !== undefined ?
        `batch.${auto!.batch[name]}` : auto!.compound[name] !== undefined ?
          `compound.${auto!.compound[name]}` : undefined);
  }

  private _mapping(): RegistrationMapping {
    const mapping: RegistrationMapping = {structure: '', compound: {},
      ...(this._target.value.peek() === 'Batches' ? {batch: {}} : {})};
    for (const [target, name] of this._columns.targets()) {
      const [prefix, field] = target.split('.');
      if (target === STRUCTURE)
        mapping.structure = name;
      else if (prefix === 'compound')
        mapping.compound[name] = field;
      else if (prefix === 'batch' && mapping.batch !== undefined)
        mapping.batch[name] = field;
    }
    return mapping;
  }

  private _problem(): string | null {
    this._columns.changed.value;
    const twice = this._columns.duplicates()[0];
    return twice !== undefined ? `Two columns of the file are mapped to ${this._labels.get(twice)}: keep one` :
      !this._columns.targets().has(STRUCTURE) ? 'Map a column to Structure' : null;
  }
}
