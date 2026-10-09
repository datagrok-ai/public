/** Project discussion and Ask SME in the app (F-51): Discuss starts or continues a thread about the
 * program, a criterion, a run or a document; Ask SME puts a question with its context to the
 * people the portfolio names. Both post into the platform chat on the program's folder. Every
 * control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divH, divV, span} from '@datagrok-libraries/u2';
import {PortfolioRow, ProgramRow} from '../generated/db';
import {askSme, askSmeOf, postToProgram} from '../discussion';
import {Anchor, describeAnchor, smeContextLines, SmeContext, smeTargetOf} from '../model/discussion';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const field = (label: string, control: HTMLElement, name: string): HTMLElement =>
  named(divH([span(label, 'pp-field-label'), control], 'pp-field'), name);

/** A post to the program's discussion, about [anchor] when given. */
export function openDiscussDialog(opts: {program: ProgramRow, anchor?: Anchor | null, onPosted?: () => void | Promise<void>}): DG.Dialog {
  installStyles();
  const text = named(Object.assign(document.createElement('textarea'), {rows: 4,
    placeholder: 'What you want to discuss: a question, a finding, a proposal'}), 'discuss-text-input');
  const errors = named(divV([], 'pp-warnings'), 'discuss-errors');
  const about = opts.anchor ? `About ${describeAnchor(opts.anchor)}` : `About ${opts.program.code} in general`;
  const dlg = ui.dialog({title: `Discuss · ${opts.program.code}`});
  dlg.root.dataset.u2Name = 'discuss-dialog';
  dlg.add(named(divV([
    named(span(about, 'pp-muted'), 'discuss-about'),
    text,
    span('Posted to the program\'s discussion: the platform chat on its workspace folder, read by whoever may open the folder.',
      'pp-muted'),
    errors,
  ], 'pp-save-run'), 'discuss-form'));
  const post = () => dlg.getButton('POST');
  dlg.addButton('POST', async () => {
    if (post()?.classList.contains('pp-busy'))
      return;
    post()?.classList.add('pp-busy', 'disabled');
    try {
      await postToProgram(opts.program, text.value, {anchor: opts.anchor ?? null});
      dlg.close();
      grok.shell.info(`${opts.program.code}: posted to the discussion`);
      await opts.onPosted?.();
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'discuss-error'));
      post()?.classList.remove('pp-busy', 'disabled');
    }
  });
  dlg.show({width: 520});
  named(post(), 'post-discussion-button');
  return dlg;
}

/** Ask SME: the question, the context it carries (read-only), and who is asked. */
export function openAskSmeDialog(opts: {portfolio: PortfolioRow, program: ProgramRow, context: SmeContext,
  onAsked?: () => void | Promise<void>}): DG.Dialog {
  installStyles();
  const config = askSmeOf(opts.portfolio);
  const target = smeTargetOf(config);
  const who = target.channel === 'group' ? `the members of ${target.groups[0]}` :
    target.channel === 'user' ? target.logins.join(', ') : `the approvers of ${opts.portfolio.code} and of the program`;
  const question = named(Object.assign(document.createElement('textarea'), {rows: 4,
    placeholder: 'Your question, in a sentence or two: the context below goes with it'}), 'sme-question-input');
  const other = named(Object.assign(document.createElement('input'), {className: 'pp-text',
    placeholder: 'Or name someone: logins, comma-separated'}), 'sme-login-input');
  const context = named(document.createElement('pre'), 'sme-context');
  context.className = 'pp-muted';
  context.style.whiteSpace = 'pre-wrap';
  context.textContent = smeContextLines(opts.context).join('\n');
  const errors = named(divV([], 'pp-warnings'), 'sme-errors');
  const dlg = ui.dialog({title: `Ask SME · ${opts.program.code}`});
  dlg.root.dataset.u2Name = 'ask-sme-dialog';
  dlg.add(named(divV([
    ...(config?.note ? [named(span(config.note, 'pp-muted'), 'sme-note')] : []),
    question,
    span('Context, sent with the question', 'pp-heading'),
    context,
    named(span(`Asks ${who}`, 'pp-muted'), 'sme-who'),
    field('Someone else', other, 'sme-login-field'),
    ...(target.notice ? [named(span(target.notice, 'pp-gap'), 'sme-external')] : []),
    span('The question goes into the program\'s discussion. The people asked are told in the platform\'s notifications ' +
      'where the server allows it; nothing is sent outside the platform.', 'pp-muted'),
    errors,
  ], 'pp-save-run'), 'ask-sme-form'));
  const ask = () => dlg.getButton('ASK');
  dlg.addButton('ASK', async () => {
    if (ask()?.classList.contains('pp-busy'))
      return;
    ask()?.classList.add('pp-busy', 'disabled');
    try {
      const logins = other.value.split(',').map((s) => s.trim()).filter((s) => s);
      const r = await askSme(opts.program.id, question.value, opts.context, logins.length ? {logins} : undefined);
      dlg.close();
      grok.shell.info(`${opts.program.code}: ${r.line}`);
      await opts.onAsked?.();
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'sme-error'));
      ask()?.classList.remove('pp-busy', 'disabled');
    }
  });
  dlg.show({width: 560});
  named(ask(), 'sme-ask-button');
  return dlg;
}
