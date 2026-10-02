import dayjs from 'dayjs';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  approvalOf, checkProposal, checkReview, checkWithdraw, describeApproval, describeProgress, needsApproval, RequestRef,
  statusAfter, validateApproval, waitingFor,
} from '../model/approvals';
import {validatePortfolio} from '../model/config';
import {describeAutoAdvance} from '../model/auto-advance';
import {documentCards} from '../model/documents';
import {ApprovalRequestRow, portalDb, PortfolioRow, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles, recordDecision, setCriterion, stagesOf} from '../portal';
import {decideOrPropose, requestsOf, reviewRequest, submitDocument, submitModelCard, waitingForMe, withdrawRequest} from '../approvals';
import {documentsOf, linkDocumentTo, PACKAGE_FILES} from '../documents';
import {loadCatalogues} from '../catalogue';
import {autoAdvance} from '../auto-advance';
import {gpc3AtHitFinding} from './auto-advance-tests';
import {withRestrictedUser} from './restricted-user';
import {_package} from '../package-test';

const request = (over: Partial<RequestRef> = {}): RequestRef => ({id: 'r1', kind: 'decision', status: 'pending',
  title: 'EXB-101: advance to Lead Generation', proposedBy: 'ann', proposedOn: '2026-10-02T09:00:00Z', source: 'person',
  needed: 1, reviews: [], programId: 'p1', portfolioId: 'pf', ...over});
const review = (by: string, verdict: 'approve' | 'reject' = 'approve') => ({by, on: '2026-10-02T10:00:00Z', verdict, rationale: 'ok'});

category('Portal: approval rules', () => {
  test('approval is optional; a portfolio names decisions and artifacts, and how many approvers', async () => {
    expect(validateApproval(undefined, 'portal').length, 0);
    expect(validateApproval({decisions: ['advance', 'terminate'], artifacts: ['document', 'model_card'], approvers: 2}, 'portal').length, 0);
    expect(validateApproval({decisions: ['sync' as any]}, 'portal')[0], 'approval.decisions: "sync" is not one of advance, recycle, hold, terminate');
    expect(validateApproval({decisions: ['advance']}, 'feed')[0].startsWith('approval.decisions needs stageAuthority "portal"'), true);
    expect(validateApproval({artifacts: ['document']}, 'feed').length, 0, 'a feed portfolio may still ask for document approval');
    expect(validateApproval({artifacts: ['memo' as any]}, 'portal')[0], 'approval.artifacts: "memo" is not one of document, model_card');
    expect(validateApproval({decisions: ['advance'], approvers: 0}, 'portal')[0], 'approval.approvers must be a whole number from 1 to 5');
    expect(validateApproval({}, 'portal')[0], 'approval names no decisions and no artifacts; leave it out to turn the workflow off');
  });

  test('only the decisions and artifacts the portfolio names need a second person', async () => {
    const a = approvalOf({decisions: ['advance'], artifacts: ['document']});
    expect(`${needsApproval(a, 'decision', 'advance')} ${needsApproval(a, 'decision', 'hold')} ${needsApproval(a, 'document')} ` +
      `${needsApproval(a, 'model_card')} ${needsApproval(null, 'decision', 'advance')}`, 'true false true false false');
    expect(approvalOf(null), null);
    expect(approvalOf({decisions: []}), null, 'an empty config is off');
    expect(describeApproval(a), 'Needs 1 approver besides the proposer: advance, document');
    expect(describeApproval(null), 'Decisions and approvals take effect at once');
  });

  test('a proposal is refused where no approval is asked, to someone who may not, and while another waits', async () => {
    const a = approvalOf({decisions: ['advance']});
    expect(checkProposal({kind: 'decision', config: a, decision: 'advance', mayPropose: true, pending: []}).length, 0);
    expect(checkProposal({kind: 'decision', config: a, decision: 'hold', mayPropose: true, pending: []})[0],
      'the portfolio does not ask for approval of hold: record it');
    expect(checkProposal({kind: 'decision', config: a, decision: 'advance', mayPropose: false, pending: []})[0],
      'Only approvers propose gate decisions');
    expect(checkProposal({kind: 'decision', config: a, decision: 'advance', mayPropose: true, pending: [request()]})[0],
      '"EXB-101: advance to Lead Generation" is waiting for approval: approve, reject or withdraw it first');
    expect(checkProposal({kind: 'decision', config: a, decision: 'advance', mayPropose: true,
      pending: [request({status: 'rejected'})]}).length, 0, 'a closed request does not block');
  });

  test('the proposer never approves their own proposal; a proposal by data any approver may', async () => {
    const ok = {reviewer: 'bob', isApprover: true, verdict: 'approve', rationale: 'data are clear'};
    expect(checkReview(request(), ok).length, 0);
    expect(checkReview(request(), {...ok, reviewer: 'ann'})[0], 'You proposed it: another approver must review it');
    expect(checkReview(request({source: 'data'}), {...ok, reviewer: 'ann'}).length, 0, 'data proposed it');
    expect(checkReview(request(), {...ok, isApprover: false})[0], 'Only approvers of the portfolio or the program review a request');
    expect(checkReview(request(), {...ok, rationale: ' '})[0], 'A review needs a rationale');
    expect(checkReview(request({reviews: [review('bob')]}), ok)[0], 'You have reviewed it already');
    expect(checkReview(request({status: 'approved'}), ok)[0], '"EXB-101: advance to Lead Generation" is approved, not waiting for approval');
    expect(checkWithdraw(request(), 'ann').length, 0);
    expect(checkWithdraw(request(), 'bob')[0], 'Only the person who proposed it withdraws it; an approver rejects it');
  });

  test('N approvers take effect together; one rejection ends it', async () => {
    expect(statusAfter([review('bob')], 1), 'approved');
    expect(statusAfter([review('bob')], 2), 'pending');
    expect(statusAfter([review('bob'), review('bob')], 2), 'pending', 'the same person twice counts once');
    expect(statusAfter([review('bob'), review('cy')], 2), 'approved');
    expect(statusAfter([review('bob'), review('cy', 'reject')], 2), 'rejected');
    expect(describeProgress({reviews: [review('bob')], needed: 2, status: 'pending'}), '1 of 2 approvals');
  });

  test('waiting for me: pending, not mine unless by data, not reviewed by me, where I approve', async () => {
    const rs = [request({id: 'a'}), request({id: 'b', proposedBy: 'me'}), request({id: 'c', source: 'data', proposedBy: 'me'}),
      request({id: 'd', reviews: [review('me')], needed: 2}), request({id: 'e', status: 'approved'}),
      request({id: 'f', portfolioId: 'other'})];
    expect(waitingFor(rs, 'me', (r) => r.portfolioId === 'pf').map((r) => r.id).join(), 'a,c');
  });

  test('advance by data in a portfolio that asks for approval is a proposal, and a submitted document waits', async () => {
    const from = {id: 's1', code: 'HF', name: 'Hit Finding', ordinal: 1};
    const to = {id: 's2', code: 'H2L', name: 'Hit-to-Lead', ordinal: 2};
    expect(describeAutoAdvance({program: 'X', from, to, reason: '', evidence: '', proposal: 'new'}),
      'X: HF → H2L by data, proposed for approval');
    expect(describeAutoAdvance({program: 'X', from, to, reason: '', evidence: '', proposal: 'pending'}),
      'X: HF → H2L by data, already proposed and waiting for approval');
    const cards = documentCards([{id: 'c1', stageId: 's2', code: 'dev', name: 'Developability', kind: 'document', required: true}],
      {current: null, next: {id: 's2', name: 'Lead'}}, [{id: 'd1', kind: 'dev', title: 'Report', status: 'submitted',
        criterionId: 'c1', path: 'x'}], {});
    expect(cards[0].missing, 'submitted: waiting for approval');
  });

  test('EX-BIOD asks for approval of advance, terminate, documents and model cards; the others do not', async () => {
    const files = await readPortfolioFiles(_package);
    const biod = files.find((c) => c.code === 'EX-BIOD')!;
    expect(JSON.stringify(biod.approval), '{"decisions":["advance","terminate"],"artifacts":["document","model_card"],"approvers":1}');
    expect(validatePortfolio(biod).join('; '), '');
    expect(files.filter((c) => c.code !== 'EX-BIOD').every((c) => !c.approval), true);
  });
});

/** A manual EX-BIOD program at [stage], removed with what it gathered afterwards. */
async function scratchProgram<T>(stageCode: string, fn: (p: ProgramRow, biod: PortfolioRow, stage: (c: string) => string) => Promise<T>):
  Promise<T> {
  const biod = await portfolioByCode('EX-BIOD');
  const stages = await stagesOf(biod.stage_model_id);
  const stage = (c: string) => stages.find((s) => s.code === c)!.id;
  const code = `TEST-AP-${Date.now() % 1e8}`;
  const [r] = await portalDb.programs.insert({portfolio_id: biod.id, code, name: code, stage_model_id: biod.stage_model_id,
    stage_id: stage(stageCode), source: 'manual', modality: 'MAB'});
  try {
    return await fn(await portalDb.programs.get(r.id), biod, stage);
  } finally {
    for (const q of await requestsOf(r.id).catch(() => []))
      await portalDb.approvalRequests.delete(q.id).catch(() => undefined);
    for (const d of await documentsOf(r.id).catch(() => []))
      await portalDb.documents.delete(d.id).catch(() => undefined);
    await portalDb.programs.delete(r.id);
  }
}

/** The request as if another person had proposed it: the proposer cannot approve their own, and
 * a second approver who can read the program row waits on GROK-21041. Withdraws [r] and files the
 * same request again under [otherUserId], with its target put back to waiting. */
async function asProposedBy(r: ApprovalRequestRow, otherUserId: string): Promise<ApprovalRequestRow> {
  await withdrawRequest(r.id);
  if (r.kind === 'document')
    await portalDb.documents.update(r.target!, {status: 'submitted'});
  if (r.kind === 'model_card')
    await portalDb.modelCards.update(r.target!, {approval: 'in_review'});
  const {id: _id, version: _v, created_on: _c, updated_on: _u, author_id: _a, ...rest} = r as any;
  const [copy] = await portalDb.approvalRequests.insert({...rest, proposed_by: otherUserId, proposed_on: dayjs(), status: 'pending'});
  return portalDb.approvalRequests.get(copy.id);
}

const refusal = async (fn: () => Promise<unknown>): Promise<string> => {
  try {
    await fn();
    return 'not refused';
  } catch (e: any) {
    return e?.message ?? String(e);
  }
};

category('Portal: approval workflow', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  test('an advance in EX-BIOD is proposed, refused to its proposer and to a non-approver, and recorded once another approves', async () => {
    await scratchProgram('TI', async (p, biod, stage) => {
      const lg = (await portalDb.gateCriterions.query({filter: `stage_id = "${stage('LG')}"`})).find((c) => c.code === 'binding')!;
      await setCriterion(p.id, lg.id, 'met');
      expect(await refusal(() => recordDecision(p.id, 'advance', null, 'binders found')),
        'EX-BIOD asks for approval of advance: propose it, and another approver approves it');
      const out = await decideOrPropose(p.id, 'advance', 'binders found');
      expect(out.recorded === undefined, true, 'nothing recorded yet');
      expect(out.proposed!.line, `${p.code}: advance to Lead Generation: proposed; waits for 1 approval`);
      expect((await portalDb.programs.get(p.id)).stage_id, stage('TI'), 'the program has not moved');
      expect(await refusal(() => decideOrPropose(p.id, 'terminate', 'stop')),
        `"${p.code}: advance to Lead Generation" is waiting for approval: approve, reject or withdraw it first`);
      expect(await refusal(() => reviewRequest(out.proposed!.request.id, 'approve', 'fine')),
        'You proposed it: another approver must review it');
      const ran = await withRestrictedUser('ppappr', async (user) => {
        // the server answers as the user; DG.User.current() stays the session's, so only the first reason is theirs
        expect((await user.asUser(() => refusal(() => reviewRequest(out.proposed!.request.id, 'approve', 'fine')))).split('; ')[0],
          'Only approvers of the portfolio or the program review a request');
        const theirs = await asProposedBy(out.proposed!.request, user.id);
        expect((await waitingForMe()).some((r) => r.id === theirs.id), true, 'waiting for admin');
        const done = await reviewRequest(theirs.id, 'approve', 'binders confirmed');
        expect(done.line, `${p.code}: advance to Lead Generation: advance recorded with 1 piece(s) of evidence`, 'the readiness frozen at the proposal');
        const after = await portalDb.approvalRequests.get(theirs.id);
        expect(`${after.status} ${!!after.stage_transition_id}`, 'approved true');
        expect((await portalDb.programs.get(p.id)).stage_id, stage('LG'));
        const t = await portalDb.stageTransitions.get(after.stage_transition_id!);
        expect(`${t.decision} ${t.source} ${t.rationale}`, 'advance approver binders found');
        expect(done.notified, 'notified', 'the proposer is told');
        return true;
      });
      if (ran === null)
        console.log('skipped the two-person half: no self-signup on this stand');
      expect(biod.code, 'EX-BIOD');
    });
  });

  test('a rejected proposal is logged and changes nothing; hold is not asked for, so it is recorded at once', async () => {
    await scratchProgram('TI', async (p, _biod, stage) => {
      const lg = (await portalDb.gateCriterions.query({filter: `stage_id = "${stage('LG')}"`})).find((c) => c.code === 'binding')!;
      await setCriterion(p.id, lg.id, 'met');
      const out = await decideOrPropose(p.id, 'advance', 'binders found');
      await withRestrictedUser('ppappr', async (user) => {
        const theirs = await asProposedBy(out.proposed!.request, user.id);
        expect(await refusal(() => reviewRequest(theirs.id, 'reject', '')), 'A review needs a rationale');
        expect((await reviewRequest(theirs.id, 'reject', 'binding is not confirmed')).line, `${p.code}: advance to Lead Generation: rejected`);
        const after = await portalDb.approvalRequests.get(theirs.id);
        expect(`${after.status} ${after.result} ${(after.reviews as any).items[0].rationale}`, 'rejected rejected binding is not confirmed');
      });
      expect((await portalDb.programs.get(p.id)).stage_id, stage('TI'), 'nothing changed');
      const held = await decideOrPropose(p.id, 'hold', 'waiting for material');
      expect(`${held.proposed === undefined} ${(await portalDb.programs.get(p.id)).status}`, 'true on_hold');
    });
  });

  test('a document in EX-BIOD is submitted, not approved at once, and approved by another approver, which meets its criterion', async () => {
    await scratchProgram('LG', async (p, _biod, stage) => {
      const dev = (await portalDb.gateCriterions.query({filter: `stage_id = "${stage('LO')}"`})).find((c) => c.code === 'developability')!;
      const path = `${PACKAGE_FILES}documents/ex-bioc/exb-201-upstream-development-report-v1.0.md`;
      const linked = await linkDocumentTo(p, {path, criterion: 'developability', title: 'Developability report', version: '1.0'});
      expect(linked.document.status, 'draft');
      expect(await refusal(() => linkDocumentTo(p, {path, status: 'approved'})),
        'EX-BIOD asks for approval of documents: submit it, and another approver approves it');
      const sub = await submitDocument(p, linked.document.id, 'ready for review');
      expect(sub.line, `${p.code}: approve Developability report · v1.0: submitted for approval`);
      expect((await portalDb.documents.get(linked.document.id)).status, 'submitted');
      expect(await refusal(() => linkDocumentTo(p, {path, title: 'Renamed'})),
        '"Developability report" is waiting for approval: approve, reject or withdraw the request first');
      await withRestrictedUser('ppappr', async (user) => {
        const theirs = await asProposedBy(sub.request, user.id);
        const done = await reviewRequest(theirs.id, 'approve', 'complete');
        expect(done.line, `${p.code}: approve Developability report · v1.0: document approved; developability met`);
      });
      expect((await portalDb.documents.get(linked.document.id)).status, 'approved');
      const status = (await portalDb.criterionStatuses.query({filter: `program_id = "${p.id}" and gate_criterion_id = "${dev.id}"`}))[0];
      expect(status.status, 'met');
    });
  });

  test('a model card is submitted and approved; reloading the catalogue keeps the approval', async () => {
    const code = `TESTAP${Date.now() % 1e6}:Model`;
    const catalogue = {models: [{code, name: 'Approval test model', card: {version: '1.0', question: 'Does it?', approval: 'draft' as const}}]};
    await loadCatalogues([{portfolio: 'EX-BIOD', catalogue}]);
    const model = (await portalDb.models.getByKey({code}))!;
    const card = (await portalDb.modelCards.getByKey({model_id: model.id, model_version: '1.0'}))!;
    try {
      const biod = await portfolioByCode('EX-BIOD');
      const smd = await portfolioByCode('EX-SMD');
      expect(await refusal(() => submitModelCard(smd, code)),
        `${code} is not in EX-SMD's catalogue`);
      const sub = await submitModelCard(biod, code, null, 'validated on the holdout set');
      expect(`${sub.line} ${(await portalDb.modelCards.get(card.id)).approval}`,
        'Approval test model 1.0: approve model card: submitted for approval in_review');
      await withRestrictedUser('ppappr', async (user) => {
        const theirs = await asProposedBy(sub.request, user.id);
        expect((await reviewRequest(theirs.id, 'approve', 'card is complete')).line, 'Approval test model 1.0: approve model card: model card approved');
      });
      const approved = await portalDb.modelCards.get(card.id);
      expect(`${approved.approval} ${!!approved.approved_on}`, 'approved true');
      const report = await loadCatalogues([{portfolio: 'EX-BIOD', catalogue}]);
      expect((await portalDb.modelCards.get(card.id)).approval, 'approved', 'the file\'s draft does not undo it');
      expect(report.includes('approval kept as approved'), true, report);
    } finally {
      for (const q of await portalDb.approvalRequests.query({filter: `target = "${card.id}"`}))
        await portalDb.approvalRequests.delete(q.id).catch(() => undefined);
      await portalDb.modelCards.delete(card.id).catch(() => undefined);
      await portalDb.models.delete(model.id).catch(() => undefined);
    }
  });

  test('advance by data in a portfolio that asks for approval proposes, waits, and moves once approved', async () => {
    const smd = await portfolioByCode('EX-SMD');
    await portalDb.portfolios.update(smd.id, {approval: {decisions: ['advance']}} as any);
    try {
      await gpc3AtHitFinding(async (p, stage) => {
        const first = await autoAdvance('EX-SMD', false);
        expect(describeAutoAdvance(first.outcomes.find((o) => o.program === 'EXS-GPC3')!), 'EXS-GPC3: HF → H2L by data, proposed for approval');
        expect((await portalDb.programs.get(p.id)).stage_id, stage('HF'), 'not moved');
        const again = await autoAdvance('EX-SMD', false);
        expect(again.outcomes.find((o) => o.program === 'EXS-GPC3')!.proposal, 'pending', 'one proposal, not two');
        const r = (await requestsOf(p.id)).find((x) => x.status === 'pending')!;
        expect(`${r.source} ${(r.evidence as any).text.includes('validated_series')}`, 'data true');
        // proposed by data, so the approver who ran it may approve
        expect((await reviewRequest(r.id, 'approve', 'data are clear')).line,
          'EXS-GPC3: advance to Hit-to-Lead: advance recorded with 2 piece(s) of evidence');
        expect((await portalDb.programs.get(p.id)).stage_id, stage('H2L'));
        const t = await portalDb.stageTransitions.get((await portalDb.approvalRequests.get(r.id)).stage_transition_id!);
        expect(`${t.source} ${t.evidence!.includes('potency')}`, 'data true');
      });
    } finally {
      await portalDb.portfolios.update(smd.id, {approval: null} as any);
      // EX-SMD asks for no approvals: every request about EXS-GPC3 is this test's
      const gpc3 = (await portalDb.programs.query({filter: `portfolio_id = "${smd.id}" and code = "EXS-GPC3"`}))[0];
      for (const q of await requestsOf(gpc3.id))
        await portalDb.approvalRequests.delete(q.id).catch(() => undefined);
    }
    expect((await portfolioByCode('EX-SMD')).approval ?? null, null, 'EX-SMD put back');
  });
});

