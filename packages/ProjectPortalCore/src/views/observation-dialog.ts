/** Record an observation (F-26): what it is about (the program, a series, campaign or study, a
 * compound), a title, the finding, and the evidence it rests on (the program's saved runs and
 * snapshots, a link). Retracting asks for a reason. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import {divV, span} from '@datagrok-libraries/u2';
import {PortfolioRow, ProgramRow} from '../generated/db';
import {addObservation, observationTargetsOf, retractObservation} from '../observations';
import {ObservationInput, ObservationSubject, ObservationView, subjectOptions} from '../model/observations';
import {linkEvidence} from '../model/evidence';
import {savedRunsOf} from '../model-runs';
import {snapshotsOf} from '../snapshots';
import {area, check, choice, choiceValue, ErrorRegion, named, openDialog, PortalDialog, text} from './u2-controls';

export interface ObservationDialogOptions {
  portfolio: PortfolioRow;
  program: ProgramRow;
  /** Preselected subject, e.g. from a series row. */
  subject?: {subject: ObservationSubject, code: string | null};
  onAdded?: (o: ObservationView) => void | Promise<void>;
}

export async function openObservationDialog(opts: ObservationDialogOptions): Promise<PortalDialog> {
  const {program} = opts;
  const options = subjectOptions(await observationTargetsOf(opts.portfolio, program.id));
  const key = (s: ObservationSubject, c: string | null) => `${s}|${c ?? ''}`;
  // evidence: the program's saved runs and snapshots, and a link
  const runs = await savedRunsOf(program.id).catch(() => []);
  const snaps = await snapshotsOf(program.id).catch(() => []);
  const picks = [
    ...runs.filter((r) => r.run.run_id).map((r) => ({kind: 'model_run' as const, target: r.run.run_id!,
      name: r.run.title ?? r.run.run_id!, label: `Run: ${r.run.title ?? r.run.run_id}`})),
    ...snaps.map((s) => ({kind: 'snapshot' as const, target: s.id, name: s.label, label: `Snapshot: ${s.label}`})),
  ];
  const errors = new ErrorRegion('observation-errors', 'observation-error');
  let added: ObservationView;
  return openDialog({title: `Observation · ${program.code}`, name: 'observation-dialog', width: 560, errors}, () => {
    const subject = choice('observation-subject-input', 'About',
      options.map((o) => ({code: key(o.subject, o.code), name: o.label})),
      opts.subject ? key(opts.subject.subject, opts.subject.code) : '');
    const title = text('observation-title-input', 'Title', {placeholder: 'One line: what was concluded'});
    const finding = area('observation-finding-input', 'Finding', {placeholder: 'The finding, and why', rows: 4});
    const boxes = picks.map((p, i) => ({p, box: check(`observation-evidence-${i}`, p.label)}));
    const link = text('observation-link-input', 'Link', {placeholder: 'A link or file path (optional)'});
    return {
      content: named(divV([
        subject.root, title.root, finding.root,
        ...(boxes.length ? [span('Evidence', 'pp-heading'),
          named(divV(boxes.map((b) => b.box.root)), 'observation-evidence')] : []),
        link.root,
        errors.root,
      ], 'pp-save-run'), 'observation-form'),
      buttons: [{label: 'RECORD', name: 'record-observation-button', primary: true,
        run: async () => {
          const [s, c] = choiceValue(subject).split('|');
          const input: ObservationInput = {subject: s as ObservationSubject, code: c || null, title: title.value.value,
            finding: finding.value.value,
            evidence: [...boxes.filter((b) => b.box.value.value)
              .map((b) => ({kind: b.p.kind, target: b.p.target, name: b.p.name})),
            ...(link.value.value.trim() ? [linkEvidence(link.value.value)] : [])]};
          added = await addObservation(program.id, input);
        },
        then: async () => {
          grok.shell.info(`${program.code}: observation recorded`);
          await opts.onAdded?.(added);
        }}],
    };
  });
}

/** Asks why, then retracts. */
export function openRetractDialog(o: ObservationView, onDone?: () => void | Promise<void>): PortalDialog {
  const errors = new ErrorRegion('retract-errors', 'retract-error');
  return openDialog({title: `Retract: ${o.title}`, name: 'retract-observation-dialog', width: 480, errors}, () => {
    const reason = area('retraction-input', 'Reason', {placeholder: 'Why it no longer holds', rows: 2});
    return {
      content: divV([span('Kept and marked retracted, with who, when and why.', 'pp-muted'), reason.root, errors.root],
        'pp-save-run'),
      buttons: [{label: 'RETRACT', name: 'retract-button', primary: true,
        run: () => retractObservation(o.id, reason.value.value), then: () => onDone?.()}],
    };
  });
}
