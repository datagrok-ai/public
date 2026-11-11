import * as DG from 'datagrok-api/dg';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Scope, signal} from '@datagrok-libraries/u2';

import {combinationHint, describeArm, findCombinations, parseArmName, therapyOf} from '../arm-name';
import {ArmsInput, armsFrame, armsOf} from '../arms-input';

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
    // the X1 of a compound's name is not a repeat count; the + after CD8 joins when another antibody follows
    expect(parseArmName('Cpd-X1 5 mg/kg')[0].drug, 'Cpd-X1');
    expect(parseArmName('Cpd-X1 5 mg/kg x3')[0].schedule, 'x3');
    expect(findCombinations(['Cpd-X1', 'Cpd-X2', 'Cpd-X1 + Cpd-X2']).length, 1);
    expect(parseArmName('anti-CD8+anti-PD-1').map((p) => p.drug).join(), 'anti-CD8,anti-PD-1');
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
    const kind = (n: string) => {
      const r = therapyOf(n);
      return r.known ? r.therapy : 'unknown';
    };
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
    expect(kind('Cpd-X kinase inhibitor 30 mg/kg p.o.'), 'unknown', 'a small molecule; how it acts, not said');
    expect(kind('Cpd-X 5 mg/kg'), 'unknown', 'a code: small molecule, and the app says it could not tell');
    expect(therapyOf('Cpd-X 5 mg/kg').therapy, 'small_molecule');
    // a combination takes its most variable drug; one drug not placed leaves it not known, with the placed one's kind
    expect(kind('cisplatin 5 mg/kg + anti-PD-1'), 'immunotherapy');
    expect(kind('Cpd-X + cisplatin'), 'unknown');
    expect(therapyOf('Cpd-X + anti-PD-1').placed, 'immunotherapy');
    // radiation is no drug: its companions give the kind; alone, it acts on the tumor
    expect(kind('cisplatin 5 mg/kg + 8 Gy'), 'small_molecule');
    expect(kind('anti-PD-1 + 8 Gy'), 'immunotherapy');
    expect(kind('Cpd-A + RT'), 'unknown');
    expect(kind('8 Gy'), 'small_molecule');
    expect(therapyOf('8 Gy').radiation, true);
    // whole words: "STING", and a "si" code with its number
    expect(kind('Cpd-X testing dose'), 'unknown');
    expect(kind('single dose 5 mg/kg'), 'unknown');
  });
});

category('Arm cards', () => {
  test('the name is the title; the kind read from it, asked in a syngeneic model; TGI or T/C set once', async () => {
    const root = Scope.root('arm cards test');
    let shown: HTMLElement | undefined;
    try {
      const host = signal<string | null>('immunocompetent');
      const metric = signal<'TGI' | 'TC'>('TGI');
      const input = root.build(() => new ArmsInput({host, metric}));
      // in the document: a link outside it does not get clicks
      document.body.append(shown = input.root);
      input.value = armsFrame(['anti-PD-1 10 mg/kg', 'Cpd-X 5 mg/kg p.o.'].map((name) =>
        ({name, therapy: 'small_molecule' as const, metric: 'TGI' as const, effect: 50,
          therapySource: 'name' as const})));
      await delay(50);
      const cards = [...input.root.querySelectorAll<HTMLElement>('.pmx-arm')];
      const grey = (c: HTMLElement) => c.querySelector<HTMLElement>('[data-pmx=understood]')!;
      const therapy = (i: number) => armsOf(input.value)[i].therapy;
      expect((cards[0].querySelector('.pmx-arm-title input') as HTMLInputElement).value, 'anti-PD-1 10 mg/kg');
      expect(cards[0].textContent!.includes('Arm 1'), false, 'no "Arm 1" over the name');
      expect(therapy(0), 'immunotherapy', 'read from the name');
      expect(grey(cards[0]).textContent!.includes('antibody or protein'), true);
      // a code in a syngeneic model: the question, through the immune system until answered
      expect(therapy(1), 'immune_agonist');
      expect(grey(cards[1]).textContent!.includes('works on the tumor'), true);
      const link = (c: HTMLElement, text: string) =>
        [...grey(c).querySelectorAll<HTMLElement>('.pmx-hint-link')].find((e) => e.textContent === text)!;
      link(cards[1], 'directly').click();
      await delay(50);
      expect(therapy(1), 'small_molecule', 'answered: directly');
      expect(armsOf(input.value)[1].therapySource, 'user');
      // the study's T/C, every arm
      metric.value = 'TC';
      await delay(50);
      expect(armsOf(input.value).map((a) => a.metric).join(), 'TC,TC');
      expect(cards[0].textContent!.includes('Expected T/C'), true);
      // in immunodeficient mice the code is a small molecule, said, not asked (the answer above stays the user's)
      host.value = 'immunodeficient';
      input.value = armsFrame([{name: 'Cpd-Y 5 mg/kg', therapy: 'immune_agonist', metric: 'TGI', effect: 50,
        therapySource: 'name'}]);
      await delay(50);
      expect(therapy(0), 'small_molecule');
      const card = input.root.querySelector<HTMLElement>('.pmx-arm')!;
      expect(grey(card).textContent!.includes('not recognized'), true);
      // a table without the source column (a script, an older plan): the therapy is the user's, kept
      input.value = DG.DataFrame.fromColumns([DG.Column.fromStrings('name', ['Cpd-Z']),
        DG.Column.fromStrings('therapy', ['cell_therapy']), DG.Column.fromStrings('metric', ['TGI']),
        DG.Column.fromFloat64Array('effect', Float64Array.from([50]))]);
      await delay(50);
      expect(therapy(0), 'cell_therapy');
    } finally {
      root.dispose();
      shown?.remove();
    }
  });

  test('the answer the user gave outlives a new model and name; the reading follows the model', async () => {
    const root = Scope.root('arm cards precedence test');
    let shown: HTMLElement | undefined;
    try {
      const host = signal<string | null>(null);
      const input = root.build(() => new ArmsInput({host}));
      document.body.append(shown = input.root);
      input.value = armsFrame(['Cpd-X 5 mg/kg', 'Cpd-Y', 'Cpd-Z + cisplatin'].map((name) =>
        ({name, therapy: 'small_molecule' as const, metric: 'TGI' as const, effect: 50,
          therapySource: 'name' as const})));
      await delay(50);
      const cards = [...input.root.querySelectorAll<HTMLElement>('.pmx-arm')];
      const grey = (c: HTMLElement) => c.querySelector<HTMLElement>('[data-pmx=understood]')!;
      const therapy = (i: number) => armsOf(input.value)[i].therapy;
      // no model yet (or the atlas loading): the cautious reading, and the question
      expect(therapy(1), 'immune_agonist', 'no host: through the immune system');
      // a combination with a drug not placed: the immune default for it, and the question
      expect(therapy(2), 'immune_agonist');
      expect(grey(cards[2]).textContent!.includes('works on the tumor'), true);
      [...grey(cards[0]).querySelectorAll<HTMLElement>('.pmx-hint-link')].find((e) => e.textContent === 'directly')!
        .click();
      await delay(50);
      host.value = 'immunodeficient';
      await delay(50);
      // the same cards, read again for the new model; the answer stays the user's
      expect(input.root.querySelectorAll('.pmx-arm')[1] === cards[1], true, 'the card is not rebuilt');
      expect(therapy(1), 'small_molecule', 'read again: immunodeficient');
      expect(therapy(0), 'small_molecule');
      host.value = 'immunocompetent';
      const title = cards[0].querySelector<HTMLInputElement>('.pmx-arm-title input')!;
      title.value = 'Cpd-X 10 mg/kg';
      title.dispatchEvent(new Event('input', {bubbles: true}));
      await delay(50);
      expect(therapy(0), 'small_molecule', 'the answer the user gave through a new model and a new name');
      expect(armsOf(input.value)[0].therapySource, 'user');
      expect(therapy(1), 'immune_agonist', 'read again: immunocompetent');
    } finally {
      root.dispose();
      shown?.remove();
    }
  });
});
