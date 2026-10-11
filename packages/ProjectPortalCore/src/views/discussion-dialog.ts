/** Project discussion and Ask SME in the app (F-51): Discuss starts or continues a thread about the
 * program, a criterion, a run or a document; Ask SME puts a question with its context to the
 * people the portfolio names. Both post into the platform chat on the program's folder. Text is
 * typed into u2's MessageInput: `@` picks a user, written as the platform's own mention markup, so
 * the Chats panel shows it as a link. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divV, span} from '@datagrok-libraries/u2';
import type {MessageInput} from '@datagrok-libraries/u2';
import {messageInput} from '@datagrok-libraries/u2/src/dg/index.js';
import {PortfolioRow, ProgramRow} from '../generated/db';
import {askSme, askSmeOf, postToProgram} from '../discussion';
import {Anchor, describeAnchor, mentionedUserIds, smeContextLines, SmeContext, smeTargetOf} from '../model/discussion';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

/** A compose box: u2's MessageInput with the platform's user mentions. [name] goes on the editable
 * part, `<name>-send` on its Send button. Ctrl+Enter (and Send) call [onSend]; a throw keeps the text.
 * In a dialog, [ownSend] false hides Send: the dialog's button sends instead. */
export function composer(opts: {name: string, placeholder: string, onSend: (markup: string) => Promise<void>,
  ownSend?: boolean, sendLabel?: string}): MessageInput {
  installStyles();
  const input = messageInput({placeholder: opts.placeholder, sendOn: 'ctrlEnter', onSend: opts.onSend});
  input.root.dataset.u2Name = `${opts.name}-composer`;
  const editor = input.root.querySelector<HTMLElement>('.u2-msg-editor');
  if (editor)
    editor.dataset.u2Name = opts.name;
  const send = input.root.querySelector<HTMLElement>('.u2-msg-send');
  if (send) {
    send.dataset.u2Name = `${opts.name}-send`;
    if (opts.sendLabel)
      send.textContent = opts.sendLabel;
  }
  if (opts.ownSend === false)
    input.root.classList.add('pp-msg-nosend');
  return input;
}

/** A post to the program's discussion, about [anchor] when given. */
export function openDiscussDialog(opts: {program: ProgramRow, anchor?: Anchor | null, onPosted?: () => void | Promise<void>}): DG.Dialog {
  installStyles();
  let submit = async (): Promise<void> => undefined;
  const text = composer({name: 'discuss-text-input', ownSend: false, onSend: () => submit(),
    placeholder: 'What you want to discuss: a question, a finding, a proposal. @ to mention someone'});
  const errors = named(divV([], 'pp-warnings'), 'discuss-errors');
  const about = opts.anchor ? `About ${describeAnchor(opts.anchor)}` : `About ${opts.program.code} in general`;
  const dlg = ui.dialog({title: `Discuss · ${opts.program.code}`});
  dlg.root.dataset.u2Name = 'discuss-dialog';
  dlg.add(named(divV([
    named(span(about, 'pp-muted'), 'discuss-about'),
    text.root,
    span('Posted to the program\'s discussion: the platform chat on its workspace folder, read by whoever may open the folder.',
      'pp-muted'),
    errors,
  ], 'pp-save-run'), 'discuss-form'));
  const post = () => dlg.getButton('POST');
  submit = async () => {
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
  };
  dlg.addButton('POST', () => void submit());
  dlg.onClose.subscribe(() => text.dispose());
  dlg.show({width: 520});
  text.focus();
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
  let submit = async (): Promise<void> => undefined;
  const question = composer({name: 'sme-question-input', ownSend: false, onSend: () => submit(),
    placeholder: 'Your question, in a sentence or two: the context below goes with it. @ to ask someone by name'});
  const context = named(document.createElement('pre'), 'sme-context');
  context.className = 'pp-muted';
  context.style.whiteSpace = 'pre-wrap';
  context.textContent = smeContextLines(opts.context).join('\n');
  const errors = named(divV([], 'pp-warnings'), 'sme-errors');
  const dlg = ui.dialog({title: `Ask SME · ${opts.program.code}`});
  dlg.root.dataset.u2Name = 'ask-sme-dialog';
  dlg.add(named(divV([
    ...(config?.note ? [named(span(config.note, 'pp-muted'), 'sme-note')] : []),
    question.root,
    span('Context, sent with the question', 'pp-heading'),
    context,
    named(span(`Asks ${who}, or the people you @mention instead`, 'pp-muted'), 'sme-who'),
    ...(target.notice ? [named(span(target.notice, 'pp-gap'), 'sme-external')] : []),
    span('The question goes into the program\'s discussion. The people asked are told in the platform\'s notifications ' +
      'where the server allows it; nothing is sent outside the platform.', 'pp-muted'),
    errors,
  ], 'pp-save-run'), 'ask-sme-form'));
  const ask = () => dlg.getButton('ASK');
  submit = async () => {
    if (ask()?.classList.contains('pp-busy'))
      return;
    ask()?.classList.add('pp-busy', 'disabled');
    try {
      const markup = question.value;
      const userIds = mentionedUserIds(markup);
      const r = await askSme(opts.program.id, markup, opts.context, userIds.length ? {userIds} : undefined);
      dlg.close();
      grok.shell.info(`${opts.program.code}: ${r.line}`);
      await opts.onAsked?.();
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'sme-error'));
      ask()?.classList.remove('pp-busy', 'disabled');
    }
  };
  dlg.addButton('ASK', () => void submit());
  dlg.onClose.subscribe(() => question.dispose());
  dlg.show({width: 560});
  question.focus();
  named(ask(), 'sme-ask-button');
  return dlg;
}
