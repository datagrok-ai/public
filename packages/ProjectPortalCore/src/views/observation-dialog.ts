/** Record an observation (F-26): what it is about (the program, a series, campaign or study, a
 * compound), a title, the finding, and the evidence it rests on (the program's saved runs and
 * snapshots, a link). Retracting asks for a reason. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divH, divV, span} from '@datagrok-libraries/u2';
import {PortfolioRow, ProgramRow} from '../generated/db';
import {addObservation, observationTargetsOf, retractObservation} from '../observations';
import {ObservationInput, ObservationSubject, ObservationView, subjectOptions} from '../model/observations';
import {linkEvidence} from '../model/evidence';
import {savedRunsOf} from '../model-runs';
import {snapshotsOf} from '../snapshots';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const field = (label: string, control: HTMLElement, name: string): HTMLElement =>
  named(divH([span(label, 'pp-field-label'), control], 'pp-field'), name);

export interface ObservationDialogOptions {
  portfolio: PortfolioRow;
  program: ProgramRow;
  /** Preselected subject, e.g. from a series row. */
  subject?: {subject: ObservationSubject, code: string | null};
  onAdded?: (o: ObservationView) => void | Promise<void>;
}

export async function openObservationDialog(opts: ObservationDialogOptions): Promise<DG.Dialog> {
  installStyles();
  const {program} = opts;
  const options = subjectOptions(await observationTargetsOf(opts.portfolio, program.id));
  const key = (s: ObservationSubject, c: string | null) => `${s}|${c ?? ''}`;
  const subject = document.createElement('select');
  subject.className = 'pp-select';
  subject.append(...options.map((o) => Object.assign(document.createElement('option'), {value: key(o.subject, o.code), textContent: o.label})));
  if (opts.subject)
    subject.value = key(opts.subject.subject, opts.subject.code);
  const title = Object.assign(document.createElement('input'), {className: 'pp-text', placeholder: 'One line: what was concluded'});
  const finding = Object.assign(document.createElement('textarea'), {rows: 4, placeholder: 'The finding, and why'});
  // evidence: the program's saved runs and snapshots, and a link
  const runs = await savedRunsOf(program.id).catch(() => []);
  const snaps = await snapshotsOf(program.id).catch(() => []);
  const picks = [
    ...runs.filter((r) => r.run.run_id).map((r) => ({kind: 'model_run' as const, target: r.run.run_id!,
      name: r.run.title ?? r.run.run_id!, label: `Run: ${r.run.title ?? r.run.run_id}`})),
    ...snaps.map((s) => ({kind: 'snapshot' as const, target: s.id, name: s.label, label: `Snapshot: ${s.label}`})),
  ];
  const boxes = picks.map((p, i) => {
    const box = named(Object.assign(document.createElement('input'), {type: 'checkbox'}), `observation-evidence-${i}`);
    return {p, box, row: divH([box, span(p.label)], 'pp-unit')};
  });
  const link = Object.assign(document.createElement('input'), {className: 'pp-text', placeholder: 'A link or file path (optional)'});
  const errors = named(divV([], 'pp-warnings'), 'observation-errors');

  const dlg = ui.dialog({title: `Observation · ${program.code}`});
  dlg.root.dataset.u2Name = 'observation-dialog';
  dlg.add(named(divV([
    field('About', subject, 'observation-subject-input'),
    field('Title', title, 'observation-title-input'),
    field('Finding', finding, 'observation-finding-input'),
    ...(boxes.length ? [span('Evidence', 'pp-heading'), named(divV(boxes.map((b) => b.row)), 'observation-evidence')] : []),
    field('Link', link, 'observation-link-input'),
    errors,
  ], 'pp-save-run'), 'observation-form'));
  const ok = () => dlg.getButton('RECORD');
  dlg.addButton('RECORD', async () => {
    if (ok()?.classList.contains('pp-busy'))
      return;
    const [s, c] = subject.value.split('|');
    const input: ObservationInput = {subject: s as ObservationSubject, code: c || null, title: title.value, finding: finding.value,
      evidence: [...boxes.filter((b) => b.box.checked).map((b) => ({kind: b.p.kind, target: b.p.target, name: b.p.name})),
        ...(link.value.trim() ? [linkEvidence(link.value)] : [])]};
    ok()?.classList.add('pp-busy', 'disabled');
    try {
      const o = await addObservation(program.id, input);
      dlg.close();
      grok.shell.info(`${program.code}: observation recorded`);
      await opts.onAdded?.(o);
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'observation-error'));
      ok()?.classList.remove('pp-busy', 'disabled');
    }
  });
  dlg.show({width: 560});
  named(ok(), 'record-observation-button');
  return dlg;
}

/** Asks why, then retracts. */
export function openRetractDialog(o: ObservationView, onDone?: () => void | Promise<void>): DG.Dialog {
  installStyles();
  const reason = Object.assign(document.createElement('textarea'), {rows: 2, placeholder: 'Why it no longer holds'});
  const errors = named(divV([], 'pp-warnings'), 'retract-errors');
  const dlg = ui.dialog({title: `Retract: ${o.title}`});
  dlg.root.dataset.u2Name = 'retract-observation-dialog';
  dlg.add(divV([span('Kept and marked retracted, with who, when and why.', 'pp-muted'),
    field('Reason', reason, 'retraction-input'), errors]));
  dlg.addButton('RETRACT', async () => {
    try {
      await retractObservation(o.id, reason.value);
      dlg.close();
      await onDone?.();
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'retract-error'));
    }
  });
  dlg.show({width: 480});
  named(dlg.getButton('RETRACT'), 'retract-button');
  return dlg;
}
