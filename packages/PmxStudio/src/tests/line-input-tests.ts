import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {Scope} from '@datagrok-libraries/u2';

import {tumorModels} from '../priors';
import {FIT_NOTES, LineInput, modelsFound, tumorsOf} from '../line-input';

category('Tumor models', () => {
  test('the atlas’s models are found by class, name, tumor type or organ', async () => {
    const models = await tumorModels();
    const lines = (cls: Parameters<typeof modelsFound>[1], q: string) =>
      modelsFound(models, cls, q).flatMap((g) => g.lines);
    expect(lines(null, '').length, Object.keys(models).length, 'everything, nothing filtered');
    expect(modelsFound(models, null, '').map((g) => g.cls).join(','), 'syngeneic,cdx,pdx', 'by class');
    expect(lines('cdx', '').join(','), 'MDA-MB-231,A549');
    expect(lines('pdx', '').join(','), 'PDX-NSCLC,PDX-PDAC');
    expect(lines(null, 'pancrea').join(','), 'Panc02,PDX-PDAC');
    expect(lines(null, 'lung').join(','), 'LLC,A549,PDX-NSCLC');
    expect(lines('cdx', 'breast').join(','), 'MDA-MB-231');
    expect(tumorsOf('MC38', models['MC38']), 'MC38 syngeneic tumors', 'the protocol names the model');
    expect(tumorsOf('PDX-PDAC', models['PDX-PDAC']).startsWith('PDX-PDAC patient-derived xenografts (PDX'), true);
    expect(tumorsOf('A549', models['A549']), 'A549 cell-line xenografts (CDX, lung)');
    expect(tumorsOf('X'), 'X tumors');
  });

  test('the picker holds a line name, as the functions take it', async () => {
    const root = Scope.root('line input test');
    const input = root.build(() => new LineInput({label: 'Tumor model'}));
    expect(input.value, 'MC38');
    input.value = '4T1';
    expect(input.value, '4T1');
    root.dispose();
  });

  test('adding your own model is the dropdown\'s foot, always in view, not the last line of the list', async () => {
    const root = Scope.root('line own test');
    const input = root.build(() => new LineInput({label: 'Tumor model'}));
    try {
      const own = input.root.querySelector<HTMLElement>('.pmx-line-own')!;
      expect(own.textContent!.startsWith('+ Add your own model'), true);
      expect([own.parentElement!.classList.contains('pmx-line-panel'), !!own.closest('.pmx-line-list')].join(),
        'true,false');
    } finally {
      root.dispose();
    }
  });

  test('the fit note: beside 4T1 and LLC only', async () => {
    const root = Scope.root('line note test');
    const input = root.build(() => new LineInput({label: 'Tumor model'}));
    const shown = () => input.root.querySelector<HTMLElement>('.pmx-line-note')!.style.display !== 'none';
    try {
      for (const [line, want] of [['4T1', true], ['LLC', true], ['MC38', false], ['Hepa1-6', false]] as const) {
        input.value = line;
        expect(shown(), want, line);
      }
      expect(Object.keys(FIT_NOTES).sort().join(), '4T1,LLC');
    } finally {
      root.dispose();
    }
  });
});
