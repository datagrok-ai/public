import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {SeqPalette, SeqPaletteBase} from './seq-palettes';

export class NucleotidesPalettes extends SeqPaletteBase {
  private static chromatogram: SeqPalette;

  public static get Chromatogram(): SeqPalette {
    if (this.chromatogram === void 0) {
      this.chromatogram = new NucleotidesPalettes({
        'A': 'green',
        'C': 'blue',
        'G': 'black', // orange ?
        'T': 'red', 'U': 'red',
        'others': 'gray',
      });
    }
    return this.chromatogram;
  }
}

export class Nucleotides {
  static readonly SemType: string = 'Nucleotides';

  static readonly SemTypeMultipleAlignment: string = 'NucleotidesMultipleAlignment';

  public static Names: StringDictionary = {
    'A': 'Adenine',
    'C': 'Cytosine',
    'G': 'Guanine',
    'T': 'Thymine',
    'U': 'Uracil',
  };
}
