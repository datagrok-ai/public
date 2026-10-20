/** The save-to-program dialog (F-14, FR-53): which portfolio, program, unit operation, run type
 * and purpose a run is saved with, filled in from where it was opened and from the model's
 * catalogue; then saveModelRun, and where the run went. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {divV, signal, span} from '@datagrok-libraries/u2';
import {RUN_OPTIONS, RUN_PURPOSES, RunPurpose} from '../model/model-runs';
import {Choice, describeSavedRun, saveRunDefaults, SaveRunForm, validateSaveRun} from '../model/save-run';
import {outsideValidity} from '../model/catalogue';
import {programsOf, saveRunChoices, sessionRuns} from '../save-run';
import {saveModelRun, SavedModelRun} from '../model-runs';
import {modelCardOf} from '../catalogue';
import {choice, choiceValue, ErrorRegion, named, openDialog, PortalDialog, text} from './u2-controls';

/** A labelled u2 choice named `<name>-input`; [blank] (when not null) is the no-choice option. [fill] replaces
 * its options and value without counting as the user's change. */
function select(name: string, label: string, options: Choice[], value: string, blank: string | null,
  onChange: () => void) {
  let filling = false;
  const input = choice(`${name}-input`, label, options, value, {...(blank !== null ? {blank} : {}),
    onChanged: () => {
      if (!filling)
        onChange();
    }});
  const fill = (items: Choice[], v: string) => {
    filling = true;
    try {
      input.setItems([...(blank !== null ? [{value: '', label: blank}] : []),
        ...items.map((o) => ({value: o.code, label: o.name}))]);
      input.value.value = items.some((o) => o.code === v) ? v : '';
    } finally {
      filling = false;
    }
  };
  return {root: input.root, value: () => choiceValue(input), set: (v: string) => input.value.value = v, fill};
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

export async function openSaveRunDialog(opts: SaveRunDialogOptions = {}): Promise<PortalDialog> {
  const offered = opts.call ? [opts.call] : sessionRuns();
  let call: DG.FuncCall | null = offered[0] ?? null;
  const choices = await saveRunChoices(call, opts.context);
  const form: SaveRunForm = saveRunDefaults(choices);
  let programs: Choice[] = form.portfolio ? await programsOf(form.portfolio).catch(() => []) : [];
  const warnings = named(divV([], 'pp-warnings'), 'validity-warnings');
  const errors = new ErrorRegion('save-run-errors', 'save-run-error');
  const valid = signal(false);
  let read = (): SaveRunForm => form;
  let saved: {r: SavedModelRun, line: string};

  const refresh = () => {
    const problems = call ? validateSaveRun(read(), programs) : ['No run of a catalogued model in this session yet'];
    errors.root.replaceChildren(...problems.map((p) => span(p, 'pp-muted')));
    valid.value = problems.length === 0;
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

  const dlg = openDialog({title: 'Save to program', name: 'save-run-dialog', width: 460, errors}, () => {
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
    const runPick = opts.call ? null : select('run', 'Run',
      offered.map((c, i) => ({code: String(i), name: describeCall(c)})), '0', null, () => void pickRun());
    const runLine = opts.call ? named(span(describeCall(opts.call), 'pp-muted'), 'save-run-run') : null;
    const portfolio = select('portfolio', 'Portfolio', choices.portfolios, form.portfolio, 'Choose…',
      () => void loadPrograms(portfolio.value(), ''));
    const program = select('program', 'Program', programs, form.program, 'Choose…', refresh);
    const unitOperation = select('unit-operation', 'Unit operation', choices.unitOperations, form.unitOperation, '—',
      refresh);
    const runType = select('run-type', 'Run type', choices.runTypes, form.runType, '—', refresh);
    const purpose = select('purpose', 'Purpose', RUN_PURPOSES.map((p) => ({code: p, name: p})), form.purpose, '—',
      refresh);
    const title = text('title-input', 'Title', {placeholder: 'Default: model · program', value: form.title});
    read = () => ({portfolio: portfolio.value(), program: program.value(), unitOperation: unitOperation.value(),
      runType: runType.value(), purpose: purpose.value() as RunPurpose | '', title: title.value.value.trim(),
      description: ''});
    // another model's run: its own unit operations, run types and defaults
    const pickRun = async () => {
      call = offered[Number(runPick!.value())] ?? null;
      const mine = ++seq;
      const next = await saveRunChoices(call, opts.context).catch(() => null);
      if (!next || mine !== seq)
        return;
      const f = saveRunDefaults(next);
      portfolio.fill(next.portfolios, f.portfolio);
      unitOperation.fill(next.unitOperations, f.unitOperation);
      runType.fill(next.runTypes, f.runType);
      purpose.set(f.purpose);
      title.value.value = f.title;
      void checkValidity();
      await loadPrograms(f.portfolio, f.program);
    };
    return {
      content: named(divV([
        ...(runPick ? [offered.length ? runPick.root :
          named(span('No run of a catalogued model in this session yet. Run one (Fed-batch titer (example), say), ' +
            'then save it from here or from the balloon that follows the run.', 'pp-muted'), 'no-session-runs')] : []),
        ...(runLine ? [runLine] : []),
        portfolio.root, program.root, unitOperation.root, runType.root, purpose.root, title.root,
        warnings, errors.root,
      ], 'pp-save-run'), 'save-run-form'),
      buttons: [{label: 'SAVE', name: 'save-run-button', primary: true, enabled: valid,
        run: async () => {
          const f = read();
          const r = await saveModelRun(call!, f.portfolio, f.program, {unitOperation: f.unitOperation || undefined,
            runType: f.runType || undefined, purpose: f.purpose || undefined, title: f.title || undefined});
          const line = describeSavedRun(f.program, call!.options[RUN_OPTIONS.stage] as string | undefined, r.folder,
            r.notes);
          saved = {r, line};
        },
        then: async () => {
          grok.shell.info(saved.line);
          await opts.onSaved?.(saved.r, saved.line);
        }}],
    };
  });
  void checkValidity();
  refresh();
  return dlg;
}
