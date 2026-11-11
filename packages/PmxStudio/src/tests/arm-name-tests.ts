import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {combinationHint, describeArm, findCombinations, parseArmName, therapyOf} from '../arm-name';

// invented names only: the curated labels the parser was measured on stay on the machine that holds them
category('Arm names', () => {
  test('drugs, dose, schedule and route of each treatment', async () => {
    const [a, b] = parseArmName('Cpd-X 5 mg/kg Q7D i.p. + anti-PD-1 10 mg/kg twice weekly');
    expect(a.drug, 'Cpd-X');
    expect(JSON.stringify(a.dose), JSON.stringify({value: 5, unit: 'mg/kg'}));
    expect(a.schedule, 'Q7D');
    expect(a.route, 'i.p.');
    expect(b.drug, 'anti-PD-1');
    expect(b.schedule, 'twice weekly');
    expect(describeArm([a]), 'Cpd-X · 5 mg/kg · Q7D · i.p.');
    // one treatment: a parenthesis is a comment; the arm's own control is not a treatment; a dose in Gy is radiation
    expect(parseArmName('Cpd-Y (a TLR7 agonist + carrier)').length, 1);
    expect(parseArmName('Cpd-Y + vehicle').length, 1);
    expect(parseArmName('Cpd-Y + isotype IgG').length, 1);
    expect(parseArmName('Cpd-Y 20 ug i.t. + 8 Gy')[1].key, 'radiation');
    expect(parseArmName('anti-PD-1 in Rag1-/- hosts')[0].drug, 'anti-PD-1');
  });

  test('an arm is the combination of two others whatever their doses and spellings', async () => {
    const found = findCombinations(['Cpd-X 5 mg/kg', 'aPD-1', 'Cpd-X 5 mg/kg + anti-PD-1 10 mg/kg']);
    expect(JSON.stringify(found), JSON.stringify([{ab: 2, a: 0, b: 1}]));
    // a single arm described at length is still that drug
    expect(findCombinations(['Cpd-X knockdown', 'anti-PD-1 antibody', 'Cpd-X + anti-PD-1']).length, 1);
    // two doses of one drug, or a drug with another not in the study, are not a combination of the arms
    expect(findCombinations(['Cpd-X 5 mg/kg', 'Cpd-X 10 mg/kg', 'Cpd-X 15 mg/kg']).length, 0);
    expect(findCombinations(['Cpd-X', 'Cpd-Z', 'Cpd-X + Cpd-W']).length, 0);
    // numbered compounds are different drugs
    expect(findCombinations(['Compound 1', 'Compound 2', 'Compound 1 + Compound 2']).length, 1);
    // a triplet: a pair of the study's and a third drug
    expect(JSON.stringify(findCombinations(['Cpd-X + Cpd-Z', 'Cpd-W', 'Cpd-X + Cpd-Z + Cpd-W'])),
      JSON.stringify([{ab: 2, a: 0, b: 1}]));
  });

  test('while a name is typed: a word on combinations only when it points there and none is found yet', async () => {
    const others = ['Cpd-X', 'anti-PD-1'];
    expect(combinationHint('Cpd-X +', others, false) !== null, true);
    expect(combinationHint('combo', others, false) !== null, true);
    expect(combinationHint('anti-PD-1 high dose', others, false) !== null, true, 'another arm\'s drug');
    expect(combinationHint('Cpd-Q 3 mg/kg', others, false), null);
    expect(combinationHint('Cpd-X + anti-PD-1', others, true), null, 'found: the checkbox says it');
    expect(combinationHint('Cpd-X +', ['Cpd-X'], false), null, 'one other arm: no pair to combine');
  });

  test('the kind of treatment read from the name; when in doubt, the one that asks for more animals', async () => {
    const kind = (n: string) => { const r = therapyOf(n); return r.known ? r.therapy : 'unknown'; };
    expect(kind('anti-PD-1 10 mg/kg i.p.'), 'immunotherapy');
    expect(kind('aPD-L1'), 'immunotherapy');
    expect(kind('rmIL-12 1 ug'), 'immunotherapy');
    expect(kind('anti-Xyz1 mAb'), 'immunotherapy', 'an antibody nobody placed: the wider spread');
    expect(kind('trastuzumab 10 mg/kg'), 'biologic');
    expect(kind('Cpd-X ADC 3 mg/kg'), 'biologic');
    expect(kind('Cpd-X STING agonist i.t.'), 'immune_agonist');
    expect(kind('CAR-T 1e6 cells i.v.'), 'cell_therapy');
    expect(kind('activated CD8+ T cells'), 'cell_therapy', 'the + of CD8+ joins nothing');
    expect(kind('Cpd-X neoantigen vaccine'), 'vaccine_peptide');
    expect(kind('siXyz8 5 nmol'), 'oligo_rna');
    expect(kind('gemcitabine 50 mg/kg i.p.'), 'small_molecule');
    expect(kind('Cpd-X kinase inhibitor 30 mg/kg p.o.'), 'unknown', 'a small molecule, but how it acts the name does not say');
    expect(kind('Cpd-X 5 mg/kg'), 'unknown', 'a code: small molecule, and the app says it could not tell');
    expect(therapyOf('Cpd-X 5 mg/kg').therapy, 'small_molecule');
    // a combination takes its most variable drug
    expect(kind('cisplatin 5 mg/kg + anti-PD-1'), 'immunotherapy');
  });
});
