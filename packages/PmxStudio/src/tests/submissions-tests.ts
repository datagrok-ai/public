import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {crc32, docx, fill} from '../docx';
import {unzip} from '../xlsx-raw';
import {PlanDocFacts, formBlocks, formItems, ntsForm, studyPlanForm} from '../ethics-documents';
import {emptyMeta} from '../study-meta';
import {CTD_TITLE, CtdStudy, commonDrug, ctdBlocks, ctdItems, ctdOfSaved, ctdRows, m4sOrder, minimalEfficacious,
  overviewRows, combinationRows, ibParagraphs,
  summaryItems} from '../ctd';
import {renamedExposure} from '../exposure';
import {ARRIVE_CREDIT, arriveItems, checklistBlocks, receivesLine} from '../arrive';
import {PREPARE_CREDIT, prepareItems} from '../prepare';
import type {SavedStudy} from '../study-file';

const facts: PlanDocFacts = {title: 'MC38 plan', meta: null, line: 'MC38',
  protocol: {'Design': 'MC38; 2 groups', 'Animals': '13 to implant', 'Randomization': 'At 100 mm³, balanced',
    'Outcome and analysis': 'Growth rates by Welch’s t', 'Multiple comparisons': 'One treated group: none needed',
    'End of study': 'Day 21', 'Humane limit': '2000 mm³', 'Measurements': 'Mon, Thu: days 0, 3, 7'},
  ethics: {design: 'We request 8 animals.', animals: 'About 10% are lost.', reduction: 'Reduction: one control.'},
  end: 21, vRand: 100, humane: 2000, diameter: '1.6', perGroup: [{name: 'Control', n: 8}, {name: 'Cpd-X', n: 8}],
  implant: 18, measurements: 7, therapies: ['small molecule']};

category('Submissions', () => {
  test('docx: a zip Word reads, with the text, the table and the fields to fill', async () => {
    expect(crc32(new TextEncoder().encode('123456789')), 0xcbf43926);
    const bytes = await docx([{heading: 'Title & <more>'}, {p: ['Plain ', fill('Grant code')]},
      {rows: [['Label', 'Value'], ['Species', 'mouse\nfemale']], head: true}, {check: true, text: 'Done'}],
    {title: 'T', landscape: true});
    const files = await unzip(bytes, () => true);
    expect([...files.keys()].join(','), '[Content_Types].xml,_rels/.rels,docProps/core.xml,word/document.xml,' +
      'word/styles.xml,word/_rels/document.xml.rels');
    for (const [name, xml] of files) {
      const parsed = new DOMParser().parseFromString(xml, 'application/xml');
      expect(parsed.getElementsByTagName('parsererror').length, 0, `${name} is well-formed XML`);
    }
    const doc = files.get('word/document.xml')!;
    expect(doc.includes('Title &amp; &lt;more&gt;'), true);
    expect(doc.includes('<w:color w:val="808080"/></w:rPr><w:t xml:space="preserve">[Grant code]'), true);
    expect(doc.includes('<w:tbl>') && doc.includes('<w:tblHeader/>') && doc.includes('<w:br/>'), true);
    expect(doc.includes('☒ ') && doc.includes('w:orient="landscape"'), true);
  });

  test('study plan: what the plan fills, what is left, the details to add', async () => {
    const items = formItems(studyPlanForm(facts));
    expect(items.length, 11);
    expect(items.filter((i) => i.done).map((i) => i.item).join(', '),
      'Inclusion and exclusion criteria, Randomisation and blinding, Outcome measures and statistics',
      'the EDA diagram: always the user; blinding off is said');
    // no details yet: the animals and what each group receives come from them
    expect(items.filter((i) => i.fix === 'details').map((i) => i.item).join(', '),
      'Experimental animals, Experimental procedures');
    // a block partly filled: what the plan filled, then what is left
    expect(items[0].says.startsWith('From the plan: study title · to complete: grant code, start date, end date'), true,
      items[0].says);
    expect(items[1].says.startsWith('From the plan: species, number · to complete: strain or genotype'), true,
      items[1].says);
    const meta = {...emptyMeta(), start: '2026-12-01', owner: 'jdoe', blinded: true, animals: {species: 'mouse',
      strain: 'C57BL/6', sex: 'female', age: '8 weeks', weight: '20 g', source: 'Charles River'},
    groups: [{group: 'Cpd-X', compound: 'Cpd-X', dose: 10, unit: 'mg/kg', route: 'PO', schedule: 'QD'}]};
    const filled = formItems(studyPlanForm({...facts, meta}));
    expect(filled.filter((i) => i.fix).length, 0);
    expect(filled[1].done, true, 'the animals, all six fields and the number');
    expect(filled[0].detail!.includes('End date: 2026-12-22'), true);
    // the plan's own answers (QA W1): a group no size reaches is said so, not left to fill; the balancing on volume
    // fills the confounders (the blinding stays the user's)
    const unsized = formItems(studyPlanForm({...facts, perGroup: [{name: 'Control', n: 8}, {name: 'Cpd-X', n: null}]}));
    expect((unsized[6].says.split('to complete:')[1] ?? '').includes('sample size'), false, unsized[6].says);
    expect(unsized[6].detail!.includes('Cpd-X: over 40'), true, unsized[6].detail);
    expect(unsized[8].done, true, unsized[8].says);
    // a group named for what it receives is said once, not "name: the same"
    const named = [{group: 'Cpd-X 10 mg/kg PO QD', compound: 'Cpd-X', dose: 10, unit: 'mg/kg', route: 'PO',
      schedule: 'QD'}, {group: 'B', compound: 'Cpd-Y', dose: 5, unit: 'mg/kg', route: 'IP', schedule: 'QW'}];
    expect(receivesLine(named)!.detail, 'Cpd-X 10 mg/kg PO QD\nB: Cpd-Y 5 mg/kg IP QW');
    // the document: one table per section, the user's fields grey in brackets
    const blocks = formBlocks({title: 'Study plan', sub: 'From the plan', sections: studyPlanForm(facts),
      credit: 'Credit'});
    expect(blocks.filter((b) => 'rows' in b).length, 3);
    const doc = (await unzip(await docx(blocks), (n) => n === 'word/document.xml')).get('word/document.xml')!;
    expect(doc.includes('[Grant code]') && doc.includes('We request 8 animals.'), true);
  });

  test('non-technical summary: Part A, labels in English or Spanish, the plan\'s Reduction', async () => {
    const en = ntsForm(facts, 'en');
    const es = ntsForm(facts, 'es');
    expect(en.map((s) => s.title).join(' | '), 'Project | Objectives and predicted benefits of the project | ' +
      'Predicted harms | Application of the Three Rs');
    expect(es[0].blocks[0].name, 'Título del proyecto');
    expect(formItems(en).length, 17);
    expect(formItems(es).map((i) => i.item).join(), formItems(en).map((i) => i.item).join(), 'short names stay');
    expect(en[3].blocks[1].fields[0][1], 'About 10% are lost. One control.');
    // the species and total once the details name the species (a mouse unless changed)
    expect(en[2].blocks[2].fields[0][1] === null, true, 'no details');
    expect(ntsForm({...facts, meta: emptyMeta()}, 'en')[2].blocks[2].fields[0][1], 'mouse: 18 in total');
  });
});

category('Submissions: CTD table', () => {
  const study: CtdStudy = {id: 'PMX-2026-014', name: 'MC38 study', compounds: ['Cpd-X'], animals: 'mouse, C57BL/6',
    line: 'MC38', endpointDay: 21, correction: '',
    groups: [{name: 'Vehicle', control: true, given: '', route: '', n: 10, nadir: -2.04, lost: 0},
      {name: 'Cpd-X', control: false, given: 'Cpd-X 10 mg/kg QD', route: 'PO', n: 9, tgi: 61.6, tgiLower: 40.2,
        tgiUpper: 80.4, tc: 38.4, p: 0.0004, regressions: {complete: 1, partial: 2}, nadir: -8.2, lost: 1}]};

  test('one row per group, the study on its first; the findings as M4S tabulates them', async () => {
    const rows = ctdRows([study]);
    expect(rows.length, 3);
    expect(rows[1].slice(0, 11).join(' | '), 'PMX-2026-014 | mouse, C57BL/6; MC38 tumour model; endpoint day 21 | ' +
      'Vehicle (control) |  | 10 | Control | — | — | — | −2.0%; 0 | No');
    expect(rows[2].slice(0, 2).join(), ',', 'the study is said once');
    expect(rows[2].slice(2, 10).join(' | '), 'Cpd-X 10 mg/kg QD | PO | 9 | 62 [40, 80] | 38 | <0.001 | ' +
      '1 CR, 2 PR | −8.2%; 1');
    expect(overviewRows([study])[1].slice(0, 3).join(' | '), 'Primary pharmacodynamics: antitumour activity, ' +
      'MC38 | mouse, C57BL/6 | PO');
    // a group without its dose: a field to fill, and the checklist sends to the study's details
    const bare = {...study, compounds: [], groups: [study.groups[0], {...study.groups[1], given: '', route: ''}]};
    expect(JSON.stringify(ctdRows([bare])[2][2]), JSON.stringify(fill('Cpd-X: dose and schedule')));
    expect(ctdItems(bare).filter((i) => i.fix === 'details').map((i) => i.item).join(', '),
      'Test article, Dose, schedule and route');
    expect(ctdItems(study).filter((i) => !i.done).map((i) => i.item).join(), 'Testing facility and location');
    const doc = (await unzip(await docx(ctdBlocks([study], {sub: 'x'}), {landscape: true}),
      (n) => n === 'word/document.xml')).get('word/document.xml')!;
    expect(doc.includes(CTD_TITLE) && doc.includes('Test article: ') && doc.includes('[Vol./Section]'), true);
  });
});

category('Submissions: CTD summary', () => {
  const g = (name: string, control = false) => ({name, control, given: '', route: '', n: 8, nadir: null, lost: null});
  const st = (id: string, animals: string, line: string): CtdStudy => ({id, name: id, compounds: ['Cpd-X'],
    animals, line, endpointDay: 21, correction: '', groups: [g('Vehicle', true), g('Cpd-X')]});

  test('the drug the picked studies share, the studies in M4S order, the ones not in said', async () => {
    expect(commonDrug(['Cpd-X, anti-PD-1', 'cpd-x', 'Cpd-X']), 'Cpd-X');
    expect(commonDrug(['Cpd-X', 'Cpd-Y']) === null, true);
    expect(m4sOrder([st('B', 'rat, SD', 'X'), st('C', 'mouse, nude', 'MC38'), st('A', 'mouse, C57BL/6', 'CT26')])
      .map((s) => s.id).join(), 'A,C,B');
    const items = summaryItems([st('PMX-1', 'mouse', 'MC38')], ['Old study']);
    expect(items.map((i) => `${i.item}:${i.status}`).join(), 'PMX-1:met,Old study:todo,Testing facility and ' +
      'location:todo');
    expect(items[0].says, 'MC38, 1 treated group, day 21');
  });

  test('a saved study keeps its rows with its effect, its ID as now, its groups renamed with it', async () => {
    const kept = st('', 'mouse', 'MC38');
    const saved = {meta: {...emptyMeta(), id: 'PMX-2026-020'}, exposure: {effect: {day: 21, tgi: {'Cpd-X': 50},
      ctd: kept}}} as unknown as SavedStudy;
    expect(ctdOfSaved(saved)!.id, 'PMX-2026-020');
    expect(ctdOfSaved({...saved, exposure: {effect: {day: 21, tgi: {}}}} as SavedStudy) === null, true);
    const renamed = renamedExposure(saved.exposure, new Map([['Cpd-X', 'Cpd-X 10']]))!;
    expect(renamed.effect!.ctd!.groups.map((x) => x.name).join(), 'Vehicle,Cpd-X 10');
    expect(Object.keys(renamed.effect!.tgi).join(), 'Cpd-X 10');
  });
});

category('Submissions: CTD exposure', () => {
  const pk = (auc: number) => ({auc, aucUnit: 'h·ng/mL', aucKind: 'AUCinf', cmax: auc / 10, cmaxUnit: 'ng/mL',
    profiles: 3, pooled: false});
  const grp = (name: string, p: number, x: ReturnType<typeof pk> | {why: string}) => ({name, control: false,
    given: `${name} QD`, route: 'PO', n: 8, tgi: 50, tc: 50, p, nadir: null, lost: null, pk: x});
  const study: CtdStudy = {id: 'PMX-1', name: 'S', compounds: ['Cpd-X'], animals: 'mouse', line: 'MC38',
    endpointDay: 21, correction: '', groups: [{name: 'Vehicle', control: true, given: '', route: '', n: 8,
      nadir: null, lost: null}, grp('Cpd-X 3 mg/kg', 0.3, pk(400)), grp('Cpd-X 10 mg/kg', 0.01, pk(1500)),
    grp('Cpd-X 30 mg/kg', 0.001, pk(5000)), grp('Cpd-X 60 mg/kg', 0.001, {why: 'NCA has no exposure at its dose'})]};

  test('AUC and Cmax per group, the lowest efficacious exposure, nothing extrapolated', async () => {
    const rows = ctdRows([study]);
    expect(rows[0][10], 'AUC; Cmax');
    expect([rows[2][10], rows[5][10]].join(' | '), '400 h·ng/mL; 40 ng/mL | —');
    expect(minimalEfficacious([study])!.group.name, 'Cpd-X 10 mg/kg', 'not significant at 3 mg/kg');
    const item = ctdItems(study).find((i) => i.item === 'Exposure (NCA)')!;
    expect([item.status, item.says].join(), 'met,3 of 4 groups with AUC and Cmax; lowest efficacious: Cpd-X 10 mg/kg',
      'a group NCA cannot give is said, not owed');
    expect(item.detail, 'Cpd-X 60 mg/kg: NCA has no exposure at its dose');
    const noNca = ctdItems({...study, groups: study.groups.map((g) => ({...g, pk: undefined}))});
    expect(noNca.length, 9, 'no NCA: no row');
    expect(noNca.map((i) => i.n).join(), '1,2,3,4,5,6,7,8,9', 'numbered without a gap');
    const none = {...study, groups: study.groups.map((g) => g.control ? g : {...g, pk: {why: 'not in NCA'}})};
    expect(ctdItems(none).find((i) => i.n === 9)!.status, 'todo', 'no group with exposure');
    const doc = (await unzip(await docx(ctdBlocks([study], {sub: 'x'}), {landscape: true}),
      (n) => n === 'word/document.xml')).get('word/document.xml')!;
    expect(doc.includes('layout of 2.6.5.3') && doc.includes('none extrapolated: PMX-1, Cpd-X 60 mg/kg') &&
      doc.includes('Minimal efficacious exposure: Cpd-X 10 mg/kg QD') && !doc.includes('[object'), true);
  });
});

category('Submissions: checklists', () => {
  const arrive = {control: 'Control', arms: ['A'], plan: 'MC38 plan', meta: null, randomized: true, endpoint: 1500,
    evaluationDay: 20, protocol: 'Slopes compared', analyzed: [{group: 'Control', n: 8}, {group: 'A', n: 8}]};
  test('ARRIVE: the Essential 10 alone, or with the Recommended Set under its own heading', async () => {
    expect(arriveItems(arrive).length, 10);
    expect(arriveItems(arrive).every((i) => i.group === undefined), true);
    const full = arriveItems(arrive, true);
    expect(full.length, 21);
    expect(full.map((i) => i.n).join(','), Array.from({length: 21}, (_, i) => i + 1).join(','));
    expect([...new Set(full.map((i) => i.group))].join('|'), 'Essential 10|Recommended Set');
    // the Recommended Set is the author's; what the study records is said
    expect(full.slice(10).every((i) => i.status === 'todo'), true);
    // the humane limits, not the study's endpoint (QA: "Humane endpoint 1000 mm³" was the endpoint)
    expect(full[15].says.startsWith('Humane endpoints: tumor over 1500 mm³ or weight loss over 20%'), true);
    expect(arriveItems({...arrive, endpoint: 1000, humane: 1800, weightLimit: 15}, true)[15].says
      .startsWith('Humane endpoints: tumor over 1800 mm³ or weight loss over 15%'), true);
    // what the groups receive: a combination with its doses and routes; a decimal comma in a name is no new line
    const given = [{group: 'T 1,5 mg/kg p.o. QD', compound: 'T', dose: 1.5, unit: 'mg/kg', route: 'p.o.',
      schedule: 'QD'}, {group: 'P 10 mg/kg i.p. BIW', compound: 'P', dose: 10, unit: 'mg/kg', route: 'i.p.',
      schedule: 'BIW'}, {group: 'T 3 mg/kg p.o. QD', compound: 'T', dose: 3, unit: 'mg/kg', route: 'p.o.',
      schedule: 'QD'}, {group: 'T + P', compound: 'T', dose: 1.5, unit: 'mg/kg', route: 'p.o.', schedule: 'QD'},
    {group: 'T + P', compound: 'P', dose: 10, unit: 'mg/kg', route: 'i.p.', schedule: 'BIW'}];
    const line = receivesLine(given)!;
    expect(line.says.endsWith(' · T 1.5 mg/kg p.o. QD + P 10 mg/kg i.p. BIW'), true, line.says);
    expect(line.detail, 'T 1,5 mg/kg p.o. QD\nP 10 mg/kg i.p. BIW\nT 3 mg/kg p.o. QD\n' +
      'T + P: T 1.5 mg/kg p.o. QD + P 10 mg/kg i.p. BIW');
    expect(full[18].says.includes('MC38 plan'), true);
    expect(arriveItems({...arrive, plan: null}, true)[18].detail === undefined, true);
  });
  test('author checklist and PREPARE as Word: one table per heading, a column for the author, the credit', async () => {
    const items = arriveItems(arrive, true);
    const blocks = checklistBlocks({title: 'ARRIVE 2.0 checklist', sub: 'From the analysis', items,
      credit: ARRIVE_CREDIT, last: ['Section, line', 'where']});
    expect(blocks.filter((b) => 'rows' in b).length, 2);
    expect(blocks.filter((b) => 'heading' in b).map((b) => (b as {heading: string}).heading).join('|'),
      'ARRIVE 2.0 checklist|Essential 10|Recommended Set');
    const e10 = checklistBlocks({title: 'E10', sub: '', items: arriveItems(arrive), credit: ARRIVE_CREDIT,
      last: ['Section, line', 'where']});
    expect(e10.filter((b) => 'rows' in b).length, 1);
    const doc = (await unzip(await docx(blocks, {landscape: true}), (n) => n === 'word/document.xml'))
      .get('word/document.xml')!;
    expect(doc.includes('Section, line') && doc.includes('[where]') && doc.includes('CC0') &&
      doc.includes('Humane endpoints: tumor over 1500 mm³') && !doc.includes('[object'), true);
    const prep = checklistBlocks({title: 'PREPARE checklist', sub: '', credit: PREPARE_CREDIT, last: ['Notes',
      'notes'], items: prepareItems({line: 'MC38', protocol: 'Slopes compared', humaneLimit: 2000, lastDay: 21,
      toImplant: 60, meta: null})});
    expect(prep.filter((b) => 'rows' in b).length, 3);
    expect(prep.filter((b) => 'rows' in b).reduce((a, b) => a + (b as {rows: unknown[]}).rows.length - 1, 0), 15);
    expect(PREPARE_CREDIT.includes('CC BY 4.0'), true);
    // the timescale from day 0 to the last day (days 0 to 21 are 22 days, never "21 days"); an arm no group size
    // powers leaves the design to settle, though the protocol says it
    const facts = {line: 'MC38', protocol: 'Slopes compared', humaneLimit: 2000, lastDay: 21, toImplant: 60,
      meta: null};
    expect(prepareItems(facts)[4].says, 'From day 0 to day 21; agree costs and who does what');
    expect(prepareItems({...facts, toImplant: null})[3].status, 'todo');
    expect(prepareItems(facts)[3].status, 'met');
  });
});

category('Submissions: combinations and IB', () => {
  const g = (name: string, given: string, tgi: number, p: number, nadir: number) =>
    ({name, control: false, given, route: 'PO', n: 8, tgi, tgiLower: tgi - 20, tgiUpper: tgi + 20, p, nadir, lost: 0});
  const study: CtdStudy = {id: 'PMX-2', name: 'combo', compounds: ['Cpd-X', 'Cpd-Y'], animals: 'mouse, BALB/c',
    line: 'CT26', endpointDay: 18, correction: '', groups: [
      {name: 'Vehicle', control: true, given: '', route: '', n: 8, nadir: -1, lost: 0},
      g('Cpd-X', 'Cpd-X 10 mg/kg QD', 40, 0.01, -3), g('Cpd-Y', 'Cpd-Y 5 mg/kg Q3D', 30, 0.04, -2),
      g('Cpd-X + Cpd-Y', 'Cpd-X 10 mg/kg QD + Cpd-Y 5 mg/kg Q3D', 75, 0.0001, -9)],
    combinations: [{ab: 3, a: 1, b: 2, step: 'Better than each drug alone'}]};

  test('S9: each combination against its two drugs alone, with where it stands', async () => {
    const rows = combinationRows([study]);
    expect(rows.length, 2);
    expect(rows[1].join(' | '), 'PMX-2 | CT26 | Cpd-X 10 mg/kg QD + Cpd-Y 5 mg/kg Q3D: TGI 75% (95% CI 55 to 95), ' +
      'p <0.001 | Cpd-X 10 mg/kg QD: TGI 40% (95% CI 20 to 60), p 0.010 | Cpd-Y 5 mg/kg Q3D: TGI 30% (95% CI 10 to ' +
      '50), p 0.040 | Better than each drug alone');
    expect(combinationRows([{...study, combinations: [{ab: 3, a: 1, b: 2}]}])[1][5], '—', 'not decided');
    expect(combinationRows([{...study, combinations: undefined}]).length, 1, 'kept before: no rows');
    const items = ctdItems(study);
    expect(items.find((i) => i.item === 'Combinations (ICH S9)')!.says, '1 against its drugs alone');
    expect(ctdItems({...study, combinations: []}).some((i) => i.item === 'Combinations (ICH S9)'), false);
  });

  test('IB: the models, each drug alone, the combination, body weight; nothing not recorded', async () => {
    const p = ibParagraphs([study]);
    expect(p.length, 3, 'no NCA: no exposure sentence');
    expect(p[0], 'The in vivo antitumour activity of Cpd-X and Cpd-Y was studied in 1 study in 1 tumour model ' +
      '(CT26) in the mouse. Against the control: Cpd-X 10 mg/kg QD PO in CT26: TGI 40% (95% CI 20 to 60), p 0.010 ' +
      'on day 18; Cpd-Y 5 mg/kg Q3D PO in CT26: TGI 30% (95% CI 10 to 50), p 0.040 on day 18.');
    expect(p[1].startsWith('In combination: Cpd-X 10 mg/kg QD + Cpd-Y 5 mg/kg Q3D in CT26: TGI 75%') &&
      p[1].endsWith('; better than each drug alone.'), true);
    expect(p[2], 'The lowest mean body-weight change from day 0 in a treated group was −9.0% (Cpd-X 10 mg/kg QD + ' +
      'Cpd-Y 5 mg/kg Q3D, CT26); 0 treated animals died or were removed for weight loss.');
    // a summary of one compound (the Overview's): it is the test article, the other only given with it
    expect(ibParagraphs([study], 'Cpd-X')[0].startsWith('The in vivo antitumour activity of Cpd-X was studied'), true);
    const article = ctdBlocks([study], {sub: 'x', drug: 'Cpd-X'})
      .find((b: any) => Array.isArray(b.p) && b.p[0]?.text === 'Test article: ') as any;
    expect(article.p[1], 'Cpd-X (given with Cpd-Y in combination)');
    const doc = (await unzip(await docx(ctdBlocks([study], {sub: 'x'}), {landscape: true}),
      (n) => n === 'word/document.xml')).get('word/document.xml')!;
    expect(doc.includes('(ICH S9)') && doc.includes('Investigator') && doc.includes('better than each drug alone') &&
      !doc.includes('[object'), true);
  });
});
