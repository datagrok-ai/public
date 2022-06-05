import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import {AminoacidsPalettes} from '../aminoacids';
import {SeqPaletteBase} from '../seq-palettes';
import {NucleotidesPalettes} from '../nucleotides';

export function _testPaletteN() {
  const cpChromatogram = NucleotidesPalettes.Chromatogram;

  // TODO: Check palettes implement SeqPalette interface
  // expect(cpChromatogram instanceof SeqPalette, true);

  expect(cpChromatogram instanceof SeqPaletteBase, true);

  expect(cpChromatogram instanceof NucleotidesPalettes, true);
}

export function _testPaletteAA() {
  const cpLest = AminoacidsPalettes.Lesk;
  const cpRasMol = AminoacidsPalettes.RasMol;
  const cpGrokGroups = AminoacidsPalettes.GrokGroups;

  // TODO: Check palettes implement SeqPalette interface
  // expect(cpLest instanceof SeqPalette, true);
  // expect(cpRasMol instanceof SeqPalette, true);
  // expect(cpGrokGroups instanceof SeqPalette, true);

  expect(cpLest instanceof SeqPaletteBase, true);
  expect(cpRasMol instanceof SeqPaletteBase, true);
  expect(cpGrokGroups instanceof SeqPaletteBase, true);

  expect(cpLest instanceof AminoacidsPalettes, true);
  expect(cpRasMol instanceof AminoacidsPalettes, true);
  expect(cpGrokGroups instanceof AminoacidsPalettes, true);
}
