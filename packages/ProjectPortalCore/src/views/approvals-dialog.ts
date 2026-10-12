/** The approval workflow in the app (F-30): one request as the program page and Waiting for me show
 * it (what is proposed, by whom, the reviews so far, and Approve or Reject with a rationale, or
 * Withdraw for the proposer), and the Waiting for me dialog over every request the user may review. */
import * as grok from 'datagrok-api/grok';
import {badge, button, divH, divV, Scope, span} from '@datagrok-libraries/u2';
import {ApprovalRow, portalDb} from '../generated/db';
import {describeRequest, reviewRequest, userNames, waitingForMe, withdrawRequest} from '../approvals';
import {reviewsOf} from '../model/approvals';
import {ROW_LIMIT} from '../rows';
import {named, openDialog, PortalDialog, text} from './u2-controls';

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
export function requestItem(r: ApprovalRow, o: RequestItemOptions): HTMLElement {
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
      const rationale = text(`${o.name}-rationale`, undefined, {placeholder: 'Rationale', inline: true});
      // one review at a time: a second click would append a second review, or apply the decision twice
      let busy = false;
      const review = (verdict: 'approve' | 'reject') => () => {
        if (busy)
          return;
        busy = true;
        approve.disabled = reject.disabled = true;
        void reviewRequest(r.id, verdict, rationale.value.value)
          .then((x) => o.onDone([x.line, ...(x.notified ? [`Proposer: ${x.notified}`] : [])].join('\n'))).catch(failed)
          .finally(() => {
            busy = false;
            approve.disabled = reject.disabled = false;
          });
      };
      const approve = named(button('Approve', review('approve'), {primary: true}), `${o.name}-approve`) as HTMLButtonElement;
      approve.title = 'Approve with the rationale: on the last approval it takes effect';
      const reject = named(button('Reject', review('reject')), `${o.name}-reject`) as HTMLButtonElement;
      reject.title = 'Reject with the rationale: it does not take effect, and the proposer is told';
      actions.push(rationale.root, approve, reject);
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
export async function openWaitingDialog(opts: WaitingOptions): Promise<PortalDialog> {
  const body = named(divV([], 'pp-units'), 'waiting-list');
  const result = named(span('', 'pp-muted'), 'waiting-result');
  // each draw owns its rationale inputs; the last one goes with the dialog
  let scope: Scope | null = null;
  const draw = async () => {
    const rows = await waitingForMe();
    const names = await userNames(rows.flatMap((r) => [r.proposed_by, ...reviewsOf(r.reviews).map((x) => x.by)]));
    const portfolios = new Map((await portalDb.portfolios.query({limit: ROW_LIMIT})).map((p) => [p.id, p]));
    scope?.dispose();
    scope = new Scope();
    Scope.runWith(scope, () => body.replaceChildren(...(rows.length ? rows.map((r, i) => {
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
    }) : [named(span('Nothing is waiting for you', 'pp-muted'), 'waiting-none')])));
  };
  await draw();
  const dlg = openDialog({title: 'Waiting for me', name: 'waiting-dialog', width: 600},
    () => ({content: divV([result, body])}));
  dlg.own(() => scope?.dispose());
  return dlg;
}
