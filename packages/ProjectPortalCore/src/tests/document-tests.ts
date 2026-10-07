import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  checkDocumentChange, compareVersions, criterionAfter, defaultKindFor, documentCards, documentEvidence, DocumentRef,
  documentVocabulary, orderDocuments, supersededBy, validateDocument, validateDocumentKinds,
} from '../model/documents';
import {describeEvidence} from '../model/evidence';
import {validatePortfolio} from '../model/config';
import {portalDb, ProgramRow} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, programContext, setCriterion, stagesOf} from '../portal';
import {
  documentCardsOf, documentsOf, folderFiles, linkDocumentTo, loadExampleDocuments, PACKAGE_FILES, setDocumentStatus,
  uploadDocument,
} from '../documents';
import {evidenceItemFor} from '../evidence';
import {PortalView} from '../views/portal-view';
import {_package} from '../package-test';

const doc = (id: string, over: Partial<DocumentRef> = {}): DocumentRef => ({id, kind: 'pd', title: id, version: '1',
  status: 'draft', criterionId: 'c-pd', path: `X:${id}.md`, ...over});

category('Portal: document rules', () => {
  test('the vocabulary: declared kinds, else one per document criterion, always with other', async () => {
    const criteria = [{code: 'pd', name: 'Process description'}, {code: 'ind', name: 'IND quality'}];
    expect(documentVocabulary(undefined, criteria).map((k) => k.code).join(), 'pd,ind,other');
    expect(documentVocabulary([{code: 'report', name: 'Report'}], criteria).map((k) => k.code).join(), 'report,other');
    expect(documentVocabulary([{code: 'other', name: 'Anything'}], []).length, 1);
    expect(defaultKindFor('pd', documentVocabulary(undefined, criteria)), 'pd');
    expect(defaultKindFor('zzz', documentVocabulary(undefined, criteria)), 'other');
    expect(validateDocumentKinds(undefined).length, 0);
    expect(validateDocumentKinds([{code: 'a', name: 'A'}, {code: 'a', name: 'B'}, {code: 'b c', name: 'X'}, {code: 'd', name: ''}]).join('|'),
      'documentKinds: "a" is declared twice|documentKinds: code "b c" may hold only letters, digits, \'_\', \'.\', \'-\'|documentKinds: d has no name');
  });

  test('a document needs a file, a title and a known kind, and answers only a document criterion', async () => {
    const v = documentVocabulary(undefined, [{code: 'pd', name: 'PD'}]);
    expect(validateDocument({path: 'X:a.md', title: 'A', kind: 'pd'}, v).length, 0);
    expect(validateDocument({path: ' ', title: '', kind: 'zz'}, v).length, 3);
    expect(validateDocument({path: 'X:a', title: 'A', kind: 'pd', criterionId: 'k'}, v, {code: 'tox', kind: 'data_condition'})[0],
      'criterion tox is a data condition; documents answer document criteria');
    expect(validateDocument({path: 'X:a', title: 'A', kind: 'pd', criterionId: 'k'}, v, null)[0],
      'the criterion does not apply to this program');
  });

  test('who may link and approve; what a status does to the criterion', async () => {
    expect(checkDocumentChange({isApprover: false, isContributor: true, status: 'draft'}).length, 0);
    expect(checkDocumentChange({isApprover: false, isContributor: true, status: 'approved'})[0], 'Only approvers mark a document approved');
    expect(checkDocumentChange({isApprover: false, isContributor: false, status: 'draft'}).length, 1);
    expect(checkDocumentChange({isApprover: true, isContributor: true, status: 'approved', programStatus: 'terminated'})[0],
      'the program is terminated');
    expect(criterionAfter(undefined, 'approved'), 'met');
    expect(criterionAfter('in_progress', 'approved'), 'met');
    expect(criterionAfter('waived', 'approved'), null);
    expect(criterionAfter('met', 'approved'), null);
    expect(criterionAfter(undefined, 'draft'), 'in_progress');
    expect(criterionAfter('not_started', 'draft'), 'in_progress');
    expect(criterionAfter('met', 'draft'), null);
    expect(criterionAfter('in_progress', 'superseded'), null);
  });

  test('a newer approved version supersedes; order; versions compare as numbers', async () => {
    const docs = [doc('a', {status: 'approved', version: '1.0'}), doc('b', {status: 'approved', version: '1.1'}),
      doc('c', {status: 'approved', kind: 'other'}), doc('d', {status: 'approved', criterionId: 'c-x'}), doc('e')];
    expect(supersededBy(docs, docs[1]).join(), 'a');
    expect(compareVersions('1.10', '1.9') > 0, true);
    expect(orderDocuments([doc('s', {status: 'superseded'}), doc('d1', {version: '0.2'}), doc('d2', {version: '0.10'}),
      doc('ap', {status: 'approved'})]).map((d) => d.id).join(), 'ap,d2,d1,s');
  });

  test('cards: document criteria of the current and the next gate, what is missing, runs and citations', async () => {
    const criteria = [
      {id: 'c-tp', stageId: 's1', code: 'tp', name: 'Transfer package', kind: 'document', required: true},
      {id: 'c-dc', stageId: 's2', code: 'tox', name: 'Tox lot', kind: 'data_condition', required: true},
      {id: 'c-pd', stageId: 's2', code: 'pd', name: 'Process description', kind: 'document', required: true},
      {id: 'c-ind', stageId: 's2', code: 'ind', name: 'IND quality', kind: 'document', required: false},
      {id: 'c-late', stageId: 's3', code: 'late', name: 'Late', kind: 'document', required: true},
    ];
    const docs = [doc('tp1', {criterionId: 'c-tp', status: 'approved'}), doc('pd1', {status: 'draft'}),
      doc('pd0', {status: 'superseded'})];
    const cards = documentCards(criteria, {current: {id: 's1', name: 'One'}, next: {id: 's2', name: 'Two'}}, docs,
      {'c-tp': 'met', 'c-pd': 'in_progress', 'c-ind': 'waived'}, {'c-pd': 'Feed 4'}, {tp1: ['advance 2026-10-02']});
    expect(cards.map((c) => `${c.code}:${c.gate}:${c.state}:${c.documents.map((d) => d.id).join('+')}:${c.missing}`).join(' | '),
      'tp:current:met:tp1: | pd:next:in_progress:pd1+pd0:only drafts: approving one meets the criterion | ind:next:waived::');
    expect(cards[0].citedBy.join(), 'advance 2026-10-02');
    expect(cards[1].run, 'Feed 4');
    expect(documentCards(criteria, {current: null, next: {id: 's3', name: 'Three'}}, [], {})[0].missing, 'no document yet');
    const e = documentEvidence(doc('tp1', {status: 'approved', title: 'Transfer package', version: '1.0'}), 'tp');
    expect(`${e.kind} ${e.target} ${e.context!.document} ${e.context!.criterion}`, 'file X:tp1.md tp1 tp');
    expect(describeEvidence(e), 'File: Transfer package · v1.0 approved');
  });

  test('documentKinds in a portfolio file is checked', async () => {
    const base: any = {code: 'X', name: 'X', setting: 'development', stageAuthority: 'portal', feeds: [],
      stageModel: {code: 'M', name: 'M', setting: 'development', stages: [{code: 'A', name: 'A', ordinal: 1}]}};
    expect(validatePortfolio({...base, documentKinds: [{code: 'r', name: 'Report'}]}).length, 0);
    expect(validatePortfolio({...base, documentKinds: [{code: 'r', name: ''}]})[0], 'documentKinds: r has no name');
  });
});

const EXAMPLE = (name: string) => `${PACKAGE_FILES}documents/ex-bioc/${name}`;

async function removeDocuments(programId: string): Promise<void> {
  for (const d of await documentsOf(programId))
    await portalDb.documents.delete(d.id).catch(() => undefined);
}

category('Portal: documents', () => {
  let exb201: ProgramRow;

  before(async () => {
    await loadAllPortfolios(_package);
    const bioc = await portfolioByCode('EX-BIOC');
    exb201 = (await portalDb.programs.query({filter: `portfolio_id = "${bioc.id}" and code = "EXB-201"`}))[0];
  });

  test('the EX-BIOC example documents: linked once, and EXB-201\'s card shows its draft process description', async () => {
    const report = await loadExampleDocuments(_package);
    expect(report.includes('not linked'), false, report);
    expect((await loadExampleDocuments(_package)).split('\n').every((l) => l.endsWith('already linked')), true);
    const ctx = await programContext(exb201.id);
    const card = (await documentCardsOf(ctx)).find((c) => c.code === 'process_description')!;
    expect(card.gate, 'next');
    // other tests and generated evidence tables may add documents to the same criterion
    expect(card.documents.map((d) => `${d.title} ${d.version} ${d.status}`).includes('Process description, Antibody Alpha 0.3 draft'), true);
    const df: DG.DataFrame = await grok.functions.call('ProjectPortalCore:programDocuments', {portfolio: 'EX-BIOC', program: 'EXB-201'});
    expect(df.rowCount >= 2, true);
    expect(df.col('kind')!.toList().includes('development_report'), true);
    const v = new PortalView();
    await v.load('EX-BIOC');
    await v.showProgram('EX-BIOC', exb201.id);
    const q = (n: string) => v.detailRoot.querySelector(`[data-u2-name="${n}"]`) as HTMLElement | null;
    expect(q('document-card-process-description')?.textContent?.includes('Process description, Antibody Alpha · v0.3 · draft'), true,
      q('document-card-process-description')?.textContent ?? 'no card');
    expect(q('other-documents')?.textContent?.includes('Upstream development report'), true);
  });

  test('a draft puts its criterion in progress, approval meets it, a newer approved version supersedes, a waiver stands', async () => {
    const bioc = await portfolioByCode('EX-BIOC');
    const st = await stagesOf(bioc.stage_model_id);
    const code = `TEST-DOC-${Date.now() % 1e9}`;
    const [row] = await portalDb.programs.insert({portfolio_id: bioc.id, code, name: code, stage_model_id: bioc.stage_model_id,
      stage_id: st.find((s) => s.code === 'D4')!.id, source: 'manual', modality: 'mAb'});
    const program = await portalDb.programs.get(row.id);
    try {
      const ctx = await programContext(program.id);
      const pd = ctx.criteria.find((c) => c.code === 'process_description')!;
      const ind = ctx.criteria.find((c) => c.code === 'ind_quality')!;
      const first = await linkDocumentTo(program, {path: EXAMPLE('exb-201-process-description-v0.3.md'), criterion: 'process_description',
        title: 'PD', version: '0.3'});
      expect(`${first.document.kind} ${first.document.status} ${first.criterion?.state}`, 'process_description draft in_progress');
      const approved = await setDocumentStatus(program, first.document.id, 'approved');
      expect(`${approved.document.status} ${approved.criterion?.state} ${!!approved.document.approved_on}`, 'approved met true');
      const second = await linkDocumentTo(program, {path: EXAMPLE('exb-201-upstream-development-report-v1.0.md'),
        criterion: pd.id, kind: 'process_description', title: 'PD', version: '1.0', status: 'approved'});
      expect(second.superseded, 1);
      expect((await portalDb.documents.get(first.document.id)).status, 'superseded');
      expect((await programContext(program.id)).recorded[pd.id], 'met');
      // a waived criterion stays waived
      await setCriterion(program.id, ind.id, 'waived');
      const w = await linkDocumentTo(program, {path: EXAMPLE('exb-202-ind-quality-v0.1.md'), criterion: 'ind_quality', status: 'approved'});
      expect(`${w.criterion} ${w.notes.join()}`, 'null ind_quality stays waived');
      // no criterion of another kind, no unknown kind, no missing file, no upload without a folder
      const refusal = (p: Promise<unknown>) => p.then(() => 'accepted', (e) => String(e?.message ?? e));
      expect(await refusal(linkDocumentTo(program, {path: EXAMPLE('exb-202-ind-quality-v0.1.md'), criterion: 'tox_lot'})),
        'criterion tox_lot is a data condition; documents answer document criteria');
      expect((await refusal(linkDocumentTo(program, {path: EXAMPLE('exb-202-ind-quality-v0.1.md'), kind: 'memo'}))).startsWith('document kind "memo"'), true);
      expect(await refusal(linkDocumentTo(program, {path: `${PACKAGE_FILES}documents/none.md`})), `no file at ${PACKAGE_FILES}documents/none.md`);
      expect((await refusal(uploadDocument(program, 'x.md', 'x'))).includes('has no workspace folder'), true);
      // the document as evidence: kind file, the record in its context
      const item = await evidenceItemFor(program, 'file', second.document.id);
      expect(`${item.kind} ${item.context?.document === second.document.id} ${item.context?.status}`, 'file true approved');
      const report: string = await grok.functions.call('ProjectPortalCore:linkDocument', {portfolio: 'EX-BIOC', program: code,
        path: EXAMPLE('exb-104-transfer-package-v1.0.md'), criterion: 'transfer_package', title: 'TP'});
      expect(report, `${code}: TP · draft linked; transfer_package in_progress`);
    } finally {
      await removeDocuments(program.id);
      for (const s of await portalDb.criterionStatuses.query({filter: `program_id = "${program.id}"`}))
        await portalDb.criterionStatuses.delete(s.id).catch(() => undefined);
      await portalDb.programs.delete(program.id);
    }
  });

  test('a file uploaded into EXB-201\'s workspace folder is linked by its path; the same name is not overwritten', async () => {
    const name = `TEST-upload-${Date.now() % 1e9}.md`;
    let id: string | null = null;
    try {
      const up = await uploadDocument(exb201, name, '# test upload');
      expect(up.folder, exb201.space_id);
      expect((await folderFiles(exb201)).some((f) => f.name === name && f.path === up.path), true);
      expect(await grok.dapi.files.readAsText(up.path), '# test upload');
      const linked = await linkDocumentTo(exb201, {path: up.path, folder: up.folder, kind: 'development_report', title: 'Test upload'});
      id = linked.document.id;
      expect(`${linked.document.folder === exb201.space_id} ${linked.criterion}`, 'true null');
      const again = await uploadDocument(exb201, name, 'x').then(() => 'overwritten', (e) => String(e?.message ?? e));
      expect(again.includes('is already in'), true, again);
    } finally {
      if (id)
        await portalDb.documents.delete(id).catch(() => undefined);
      await grok.dapi.spaces.id(exb201.space_id!).files.delete(name).catch(() => undefined);
    }
  });
});
