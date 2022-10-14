import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import {_testPaletteN, _testPaletteAA} from '@datagrok-libraries/bio/src/tests/palettes-tests';

category('Palettes', () => {
  test('testPaletteN', async () => { await _testPaletteN(); });
  test('testPaletteAA', async () => { await _testPaletteAA(); });

  test('testPalettePtMe', async () => {
    const colorMeNle = bio.AminoacidsPalettes.GrokGroups.get('MeNle');
    const colorMeA = bio.AminoacidsPalettes.GrokGroups.get('MeA');
    const colorMeG = bio.AminoacidsPalettes.GrokGroups.get('MeG');
    const colorMeF = bio.AminoacidsPalettes.GrokGroups.get('MeF');

    const colorL = bio.AminoacidsPalettes.GrokGroups.get('L');
    const colorA = bio.AminoacidsPalettes.GrokGroups.get('A');
    const colorG = bio.AminoacidsPalettes.GrokGroups.get('G');
    const colorF = bio.AminoacidsPalettes.GrokGroups.get('F');

    expect(colorMeNle, colorL);
    expect(colorMeA, colorA);
    expect(colorMeG, colorG);
    expect(colorMeF, colorF);
  });
});
