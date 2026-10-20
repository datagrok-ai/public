/** Project discussion and Ask SME in the app (F-51): Discuss starts or continues a thread about the
 * program, a criterion, a run or a document; Ask SME puts a question with its context to the
 * people the portfolio names. Both post into the platform chat on the program's folder. Text is
 * typed into u2's MessageInput: `@` picks a user, written as the platform's own mention markup, so
 * the Chats panel shows it as a link. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import {divV, span} from '@datagrok-libraries/u2';
import type {MessageInput} from '@datagrok-libraries/u2';
import {composer, ErrorRegion, named, openDialog, PortalDialog} from './u2-controls';
import {PortfolioRow, ProgramRow} from '../generated/db';
import {askSme, askSmeOf, postToProgram} from '../discussion';
import {Anchor, describeAnchor, mentionedUserIds, smeContextLines, SmeContext, smeTargetOf} from '../model/discussion';

/** A post to the program's discussion, about [anchor] when given. Ctrl+Enter in the box posts as POST does. */
export function openDiscussDialog(opts: {program: ProgramRow, anchor?: Anchor | null,
  onPosted?: () => void | Promise<void>}): PortalDialog {
  const about = opts.anchor ? `About ${describeAnchor(opts.anchor)}` : `About ${opts.program.code} in general`;
  const errors = new ErrorRegion('discuss-errors', 'discuss-error');
  let text!: MessageInput;
  const dlg = openDialog({title: `Discuss · ${opts.program.code}`, name: 'discuss-dialog', width: 520, errors}, () => {
    text = composer({name: 'discuss-text-input', ownSend: false,
      onSend: () => sent(dlg.press('post-discussion-button')),
      placeholder: 'What you want to discuss: a question, a finding, a proposal. @ to mention someone'});
    return {
      content: named(divV([
        named(span(about, 'pp-muted'), 'discuss-about'),
        text.root,
        span('Posted to the program\'s discussion: the platform chat on its workspace folder, read by whoever ' +
          'may open the folder.', 'pp-muted'),
        errors.root,
      ], 'pp-save-run'), 'discuss-form'),
      buttons: [{label: 'POST', name: 'post-discussion-button', primary: true,
        run: () => postToProgram(opts.program, text.value, {anchor: opts.anchor ?? null}),
        then: async () => {
          grok.shell.info(`${opts.program.code}: posted to the discussion`);
          await opts.onPosted?.();
        }}],
    };
  });
  text.focus();
  return dlg;
}

/** Ask SME: the question, the context it carries (read-only), and who is asked. */
export function openAskSmeDialog(opts: {portfolio: PortfolioRow, program: ProgramRow, context: SmeContext,
  onAsked?: () => void | Promise<void>}): PortalDialog {
  const config = askSmeOf(opts.portfolio);
  const target = smeTargetOf(config);
  const who = target.channel === 'group' ? `the members of ${target.groups[0]}` :
    target.channel === 'user' ? target.logins.join(', ') : `the approvers of ${opts.portfolio.code} and of the program`;
  const context = named(document.createElement('pre'), 'sme-context');
  context.className = 'pp-muted';
  context.style.whiteSpace = 'pre-wrap';
  context.textContent = smeContextLines(opts.context).join('\n');
  const errors = new ErrorRegion('sme-errors', 'sme-error');
  let question!: MessageInput;
  let line = '';
  const dlg = openDialog({title: `Ask SME · ${opts.program.code}`, name: 'ask-sme-dialog', width: 560, errors}, () => {
    question = composer({name: 'sme-question-input', ownSend: false, onSend: () => sent(dlg.press('sme-ask-button')),
      placeholder: 'Your question, in a sentence or two: the context below goes with it. @ to ask someone by name'});
    return {
      content: named(divV([
        ...(config?.note ? [named(span(config.note, 'pp-muted'), 'sme-note')] : []),
        question.root,
        span('Context, sent with the question', 'pp-heading'),
        context,
        named(span(`Asks ${who}, or the people you @mention instead`, 'pp-muted'), 'sme-who'),
        ...(target.notice ? [named(span(target.notice, 'pp-gap'), 'sme-external')] : []),
        span('The question goes into the program\'s discussion. The people asked are told in the platform\'s ' +
          'notifications where the server allows it; nothing is sent outside the platform.', 'pp-muted'),
        errors.root,
      ], 'pp-save-run'), 'ask-sme-form'),
      buttons: [{label: 'ASK', name: 'sme-ask-button', primary: true,
        run: async () => {
          const markup = question.value;
          const userIds = mentionedUserIds(markup);
          line = (await askSme(opts.program.id, markup, opts.context, userIds.length ? {userIds} : undefined)).line;
        },
        then: async () => {
          grok.shell.info(`${opts.program.code}: ${line}`);
          await opts.onAsked?.();
        }}],
    };
  });
  question.focus();
  return dlg;
}

/** Ctrl+Enter in a dialog's box: a send that did not go through throws, so the box keeps the text. */
async function sent(pressed: Promise<boolean>): Promise<void> {
  if (!await pressed)
    throw new Error('not sent');
}
