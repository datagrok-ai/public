/** The save-to-program dialog (F-14, FR-53): which portfolio, program, unit operation, run type
 * and purpose a run is saved with, filled in from where it was opened and from the model's
 * catalogue; then saveModelRun, and where the run went. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';
import {divV, divH, span} from '@datagrok-libraries/u2';
import {RUN_OPTIONS, RUN_PURPOSES, RunPurpose} from '../model/model-runs';
import {Choice, describeSavedRun, saveRunDefaults, SaveRunForm, validateSaveRun} from '../model/save-run';
import {outsideValidity} from '../model/catalogue';
import {programsOf, saveRunChoices, sessionRuns} from '../save-run';
import {saveModelRun, SavedModelRun} from '../model-runs';
import {modelCardOf} from '../catalogue';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

/** A labelled native select: the BDD `selects ... in` gesture and the keyboard both work on it. */
function select(name: string, label: string, options: Choice[], value: string, blank: string | null) {
  const el = document.createElement('select');
  el.className = 'pp-select';
  const fill = (items: Choice[], v: string) => {
    el.replaceChildren(...[...(blank !== null ? [{code: '', name: blank}] : []), ...items].map((o) =>
      Object.assign(document.createElement('option'), {value: o.code, textContent: o.name})));
    el.value = items.some((o) => o.code === v) ? v : '';
  };
  fill(options, value);
  const root = named(divH([span(label, 'pp-field-label'), el], 'pp-field'), `${name}-input`);
  return {root, el, fill};
}

const describeCall = (c: DG.FuncCall): string => {
  const when = c.started ? dayjs(c.started as any).format('HH:mm:ss') : '';
  const inputs = c.func.inputs.slice(0, 3).map((p) => `${p.caption || p.name} ${c.inputs[p.name]}`).join(', ');
  return `${c.func.friendlyName ?? c.func.name}${when ? ` · ${when}` : ''}${inputs ? ` · ${inputs}` : ''}`;
};

export interface SaveRunDialogOptions {
  /** The run to save; without one the dialog offers this session's runs of catalogued models. */
  call?: DG.FuncCall | null;
  /** A program page opens the dialog on its own portfolio and program. */
  context?: {portfolio?: string | null, program?: string | null} | null;
  onSaved?: (r: SavedModelRun, line: string) => void | Promise<void>;
}

export async function openSaveRunDialog(opts: SaveRunDialogOptions = {}): Promise<DG.Dialog> {
  installStyles();
  const offered = opts.call ? [opts.call] : sessionRuns();
  let call: DG.FuncCall | null = offered[0] ?? null;
  const choices = await saveRunChoices(call, opts.context);
  const form: SaveRunForm = saveRunDefaults(choices);
  let programs: Choice[] = form.portfolio ? await programsOf(form.portfolio).catch(() => []) : [];

  const runPick = opts.call ? null : select('run', 'Run',
    offered.map((c, i) => ({code: String(i), name: describeCall(c)})), '0', null);
  const runLine = opts.call ? named(span(describeCall(opts.call), 'pp-muted'), 'save-run-run') : null;
  const portfolio = select('portfolio', 'Portfolio', choices.portfolios, form.portfolio, 'Choose…');
  const program = select('program', 'Program', programs, form.program, 'Choose…');
  const unitOperation = select('unit-operation', 'Unit operation', choices.unitOperations, form.unitOperation, '—');
  const runType = select('run-type', 'Run type', choices.runTypes, form.runType, '—');
  const purpose = select('purpose', 'Purpose', RUN_PURPOSES.map((p) => ({code: p, name: p})), form.purpose, '—');
  const title = document.createElement('input');
  title.className = 'pp-text';
  title.placeholder = 'Default: model · program';
  title.value = form.title;
  const titleRow = named(divH([span('Title', 'pp-field-label'), title], 'pp-field'), 'title-input');
  const warnings = named(divV([], 'pp-warnings'), 'validity-warnings');
  const errors = named(divV([], 'pp-warnings'), 'save-run-errors');

  const dlg = ui.dialog({title: 'Save to program'});
  dlg.root.dataset.u2Name = 'save-run-dialog';
  const save = () => dlg.getButton('SAVE');

  const read = (): SaveRunForm => ({portfolio: portfolio.el.value, program: program.el.value,
    unitOperation: unitOperation.el.value, runType: runType.el.value, purpose: purpose.el.value as RunPurpose | '',
    title: title.value.trim(), description: ''});
  const refresh = () => {
    const problems = call ? validateSaveRun(read(), programs) : ['No run of a catalogued model in this session yet'];
    errors.replaceChildren(...problems.map((p) => span(p, 'pp-muted')));
    save()?.classList.toggle('disabled', problems.length > 0);
  };
  /** Inputs outside the card's validity range: a warning, not a refusal. */
  const checkValidity = async () => {
    warnings.replaceChildren();
    const card = call ? await modelCardOf(call.func.nqName).catch(() => null) : null;
    if (!card || !call)
      return;
    const values: Record<string, unknown> = {};
    for (const p of card.inputs)
      values[p.name] = call.inputs[p.name];
    warnings.replaceChildren(...outsideValidity(card.inputs, values).map((w) => span(w, 'pp-reason')));
  };

  // a later change wins: a slow read for an earlier pick does not overwrite a newer one
  let seq = 0;
  const loadPrograms = async (code: string, keep: string) => {
    const mine = ++seq;
    const list = code ? await programsOf(code).catch(() => []) : [];
    if (mine !== seq)
      return;
    programs = list;
    program.fill(programs, keep);
    refresh();
  };
  portfolio.el.onchange = () => void loadPrograms(portfolio.el.value, '');
  for (const s of [program, unitOperation, runType, purpose])
    s.el.onchange = refresh;
  if (runPick) {
    // another model's run: its own unit operations, run types and defaults
    runPick.el.onchange = async () => {
      call = offered[Number(runPick.el.value)] ?? null;
      const mine = ++seq;
      const next = await saveRunChoices(call, opts.context).catch(() => null);
      if (!next || mine !== seq)
        return;
      const f = saveRunDefaults(next);
      portfolio.fill(next.portfolios, f.portfolio);
      unitOperation.fill(next.unitOperations, f.unitOperation);
      runType.fill(next.runTypes, f.runType);
      purpose.el.value = f.purpose;
      title.value = f.title;
      void checkValidity();
      await loadPrograms(f.portfolio, f.program);
    };
  }

  dlg.add(named(divV([
    ...(runPick ? [offered.length ? runPick.root :
      named(span('No run of a catalogued model in this session yet. Run one (Fed-batch titer (example), say), then ' +
        'save it from here or from the balloon that follows the run.', 'pp-muted'), 'no-session-runs')] : []),
    ...(runLine ? [runLine] : []),
    portfolio.root, program.root, unitOperation.root, runType.root, purpose.root, titleRow,
    warnings, errors,
  ], 'pp-save-run'), 'save-run-form'));

  dlg.addButton('SAVE', async () => {
    const f = read();
    if (!call || validateSaveRun(f, programs).length > 0 || save()?.classList.contains('pp-busy'))
      return;
    save()?.classList.add('pp-busy', 'disabled');
    try {
      const r = await saveModelRun(call, f.portfolio, f.program, {unitOperation: f.unitOperation || undefined,
        runType: f.runType || undefined, purpose: f.purpose || undefined, title: f.title || undefined});
      const line = describeSavedRun(f.program, call.options[RUN_OPTIONS.stage] as string | undefined, r.folder, r.notes);
      dlg.close();
      grok.shell.info(line);
      await opts.onSaved?.(r, line);
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'save-run-error'));
      save()?.classList.remove('pp-busy', 'disabled');
    }
  });
  dlg.show({width: 460});
  named(save(), 'save-run-button');
  void checkValidity();
  refresh();
  return dlg;
}
