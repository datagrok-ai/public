/** The approval workflow in the app (F-30): one request as the program page and Waiting for me show
 * it (what is proposed, by whom, the reviews so far, and Approve or Reject with a rationale, or
 * Withdraw for the proposer), and the Waiting for me dialog over every request the user may review. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {badge, button, divH, divV, span} from '@datagrok-libraries/u2';
import {ApprovalRequestRow, portalDb} from '../generated/db';
import {describeRequest, reviewRequest, userNames, waitingForMe, withdrawRequest} from '../approvals';
import {reviewsOf} from '../model/approvals';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

export interface RequestItemOptions {
  names: Map<string, string>;
  /** The user is an approver here and did not propose it: Approve and Reject. */
  canReview: boolean;
  /** The user proposed it: Withdraw. */
  canWithdraw: boolean;
  /** After a review or a withdrawal, with the line to show. */
  onDone: (line: string) => void | Promise<void>;
  /** The automation name of the item; its controls take it as a prefix. */
  name: string;
}

/** One request: its line, the reviews so far and, while it waits, what the user may do. */
export function requestItem(r: ApprovalRequestRow, o: RequestItemOptions): HTMLElement {
  const reviews = reviewsOf(r.reviews).map((x) =>
    span(`${o.names.get(x.by) ?? 'someone'} ${x.verdict === 'approve' ? 'approved' : 'rejected'}: ${x.rationale}`, 'pp-gap'));
  const lines: HTMLElement[] = [
    divH([badge(r.status, {variant: r.status === 'approved' ? 'success' : r.status === 'pending' ? 'warning' : 'error'}),
      named(span(describeRequest(r, o.names)), `${o.name}-line`)], 'pp-unit'),
    ...(r.rationale ? [span(`Rationale: ${r.rationale}`, 'pp-muted')] : []),
    ...reviews,
    ...(r.status !== 'pending' && r.result ? [span(r.result, 'pp-muted')] : []),
  ];
  if (r.status === 'pending' && (o.canReview || o.canWithdraw)) {
    const failed = (e: any) => grok.shell.error(e?.message ?? String(e));
    const actions: HTMLElement[] = [];
    if (o.canReview) {
      const rationale = named(document.createElement('input'), `${o.name}-rationale`);
      rationale.placeholder = 'Rationale';
      const review = (verdict: 'approve' | 'reject') => () => void reviewRequest(r.id, verdict, rationale.value)
        .then((x) => o.onDone([x.line, ...(x.notified ? [`Proposer: ${x.notified}`] : [])].join('\n'))).catch(failed);
      const approve = named(button('Approve', review('approve'), {primary: true}), `${o.name}-approve`);
      approve.title = 'Approve with the rationale: on the last approval it takes effect';
      const reject = named(button('Reject', review('reject')), `${o.name}-reject`);
      reject.title = 'Reject with the rationale: it does not take effect, and the proposer is told';
      actions.push(rationale, approve, reject);
    }
    if (o.canWithdraw) {
      const withdraw = named(button('Withdraw', () => void withdrawRequest(r.id)
        .then((x) => o.onDone(`${x.title}: withdrawn`)).catch(failed)), `${o.name}-withdraw`);
      withdraw.title = 'Take the proposal back; nothing changes';
      actions.push(withdraw);
    }
    lines.push(named(divH(actions, 'pp-actions pp-role-form'), `${o.name}-actions`));
  }
  return named(divV(lines, 'pp-attention-item'), o.name);
}

export interface WaitingOptions {
  /** Opens the program a request is about, on the portal page. */
  openProgram: (portfolioCode: string, programId: string) => Promise<void>;
}

/** Waiting for me: every pending request the current user may review, oldest first. */
export async function openWaitingDialog(opts: WaitingOptions): Promise<DG.Dialog> {
  installStyles();
  const body = named(divV([], 'pp-units'), 'waiting-list');
  const result = named(span('', 'pp-muted'), 'waiting-result');
  const dlg = ui.dialog({title: 'Waiting for me'});
  const draw = async () => {
    const rows = await waitingForMe();
    const names = await userNames(rows.flatMap((r) => [r.proposed_by, ...reviewsOf(r.reviews).map((x) => x.by)]));
    const portfolios = new Map((await portalDb.portfolios.query({})).map((p) => [p.id, p]));
    body.replaceChildren(...(rows.length ? rows.map((r, i) => {
      const item = requestItem(r, {names, canReview: true, canWithdraw: false, name: `waiting-${i}`,
        onDone: async (line) => {
          result.textContent = line;
          await draw();
        }});
      if (r.program_id) {
        const open = named(span(`Open in ${portfolios.get(r.portfolio_id)?.name ?? 'the portal'}`, 'pp-link'), `waiting-${i}-open`);
        open.onclick = () => {
          dlg.close();
          void opts.openProgram(portfolios.get(r.portfolio_id)?.code ?? '', r.program_id!);
        };
        item.append(open);
      }
      return item;
    }) : [named(span('Nothing is waiting for you', 'pp-muted'), 'waiting-none')]));
  };
  await draw();
  dlg.add(divV([result, body]));
  dlg.root.dataset.u2Name = 'waiting-dialog';
  dlg.show({width: 600});
  return dlg;
}
