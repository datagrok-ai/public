import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import {_testPaletteN, _testPaletteAA} from '@datagrok-libraries/bio/src/tests/palettes-tests';
import { AminoacidsPalettes } from '@datagrok-libraries/bio';

category('Palettes', () => {
  test('testPaletteN', async () => { await _testPaletteN(); });
  test('testPaletteAA', async () => { await _testPaletteAA(); });

  test('testPalettePtMe', async () => {
    const colorMeNle = AminoacidsPalettes.GrokGroups.get('MeNle');
    const colorMeA = AminoacidsPalettes.GrokGroups.get('MeA');
    const colorMeG = AminoacidsPalettes.GrokGroups.get('MeG');
    const colorMeF = AminoacidsPalettes.GrokGroups.get('MeF');

    const colorL = AminoacidsPalettes.GrokGroups.get('L');
    const colorA = AminoacidsPalettes.GrokGroups.get('A');
    const colorG = AminoacidsPalettes.GrokGroups.get('G');
    const colorF = AminoacidsPalettes.GrokGroups.get('F');

    expect(colorMeNle, colorL);
    expect(colorMeA, colorA);
    expect(colorMeG, colorG);
    expect(colorMeF, colorF);
  });
});
