import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {SeqPalette, SeqPaletteBase} from './seq-palettes';

export class AminoacidsPalettes extends SeqPaletteBase {
  private static lesk: SeqPalette;

  public static get Lesk(): SeqPalette {
    if (this.lesk === void 0) {
      this.lesk = this.makePalette([
        [['G', 'A', 'S', 'T'], 'orange'],
        [['C', 'V', 'I', 'L', 'P', 'F', 'Y', 'M', 'W'], 'all_green'],
        [['N', 'Q', 'H'], 'magenta'],
        [['D', 'E'], 'red'],
        [['K', 'R'], 'all_blue'],
      ], false, AminoacidsPalettes);
    }
    return this.lesk;
  }

  private static grokGroups: SeqPalette;

  public static get GrokGroups(): SeqPalette {
    if (this.grokGroups === void 0) {
      this.grokGroups = this.makePalette([
        [['C', 'U'], 'yellow'],
        [['G', 'P'], 'red'],
        [['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'], 'all_green'],
        [['R', 'H', 'K'], 'light_blue'],
        [['D', 'E'], 'dark_blue'],
        [['S', 'T', 'N', 'Q'], 'orange'],
      ], false, AminoacidsPalettes);
    }
    return this.grokGroups;
  }

  private static rasMol: SeqPalette;

  public static get RasMol(): SeqPalette {
    if (this.rasMol === void 0) {
      this.rasMol = new AminoacidsPalettes({
        // http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
        'D': '#E60A0A', // asp, aspartic acid, asp
        'E': '#E60A0A', // glu, glutamic acid
        'C': '#E6E600', // cys, cysteine
        'M': '#E6E600', // met, methionine
        'K': '#145AFF', // lys, lysine
        'R': '#145AFF', // arg, arginine
        'S': '#FA9600', // ser, serine
        'T': '#FA9600', // thr, threonine
        'F': '#3232AA', // phe, phenylalanine
        'Y': '#3232AA', // tyr, tyrosine
        'N': '#00DCDC', // asn, asparagine
        'Q': '#00DCDC', // gln, glutamine
        'G': '#EBEBEB', // gly, glycine
        'L': '#0F820F', // leu, leucine
        'V': '#0F820F', // val, valine
        'I': '#0F820F', // ile, isoleucine
        'A': '#C8C8C8', // ala, alanine
        'W': '#B45AB4', // trp, tryptophan
        'H': '#8282D2', // his, histidine
        'P': '#DC9682', // pro, proline
        'others': '#BEA06E',
      });
    }
    return this.rasMol;
  }
}

export class Aminoacids {
  static readonly SemType: string = 'Aminoacids';

  static readonly SemTypeMultipleAlignment: string = 'AminoacidsMultipleAlignment';

  static undefinedColor = 'rgb(100,100,100)';

  public static Names: StringDictionary = {
    'G': 'Glycine',
    'L': 'Leucine',
    'Y': 'Tyrosine',
    'S': 'Serine',
    'E': 'Glutamic acid',
    'Q': 'Glutamine',
    'D': 'Aspartic acid',
    'N': 'Asparagine',
    'F': 'Phenylalanine',
    'A': 'Alanine',
    'K': 'Lysine',
    'R': 'Arginine',
    'H': 'Histidine',
    'C': 'Cysteine',
    'V': 'Valine',
    'P': 'Proline',
    'W': 'Tryptophan',
    'I': 'Isoleucine',
    'M': 'Methionine',
    'T': 'Threonine',
  };

  public static AASmiles: StringDictionary = {
    'G': 'NCC(=O)O',
    'L': 'N[C@H](CC(C)C)C(=O)O',
    'Y': 'NC(CC1=CC=C(O)C=C1)C(=O)O',
    'S': 'NC(CO)C(=O)O',
    'E': 'N[C@@H](CCC(O)=O)C(=O)O',
    'Q': 'N[C@@H](CCC(N)=O)C(=O)O',
    'D': 'N[C@@H](CC(O)=O)C(=O)O',
    'N': 'N[C@@H](CC(N)=O)C(=O)O',
    'F': 'NC(CC1=CC=CC=C1)C(=O)O',
    'A': 'N[C@H](C)C(=O)O',
    'K': 'NC(CCCCN)C(=O)O',
    'R': 'N[C@H](CCCNC(=N)C)C(=O)O',
    'H': 'NC(CC1=CN=C[N]1)C(=O)O',
    'C': 'N[C@@H](CS)C(=O)O',
    'V': 'NC(C(C)C)C(=O)O',
    'P': 'N(CCC1)C1C(=O)O',
    'W': 'N[C@@H](Cc1c2ccccc2n([H])c1)C(=O)O',
    'I': 'N[C@H]([C@H](C)CC)C(=O)O',
    'M': 'NC(CCSC)C(=O)O',
    'T': 'NC(C(O)C)C(=O)O',
  };

  public static AASmilesTruncated: StringDictionary = {
    'G': '*C*',
    'L': 'CC(C)C[C@H](*)*',
    'Y': 'C1=CC(=CC=C1CC(*)*)O',
    'S': 'OCC(*)C*',
    'E': '*[C@@H](CCC(O)=O)*',
    'Q': '*N[C@@H](CCC(N)=O)*',
    'D': '*[C@@H](CC(O)=O)*',
    'N': '*[C@@H](CC(N)=O)*',
    'F': 'C1=CC=C(C=C1)CC(*)*',
    'A': 'C[C@H](*)*',
    'K': 'C(CCN)CC(*)*',
    'R': '*[C@H](CCCNC(=N)C)*',
    'H': 'C1=C(NC=N1)CC(*)*',
    'C': 'C([C@@H](*)*)S',
    'V': 'CC(C)C(*)*',
    'P': 'C1CCN(*)C1*',
    'W': '*[C@@H](Cc1c2ccccc2n([H])c1)*',
    'I': 'CC[C@H](C)[C@H](*)*',
    'M': 'CSCCC(*)*',
    'T': 'CC(O)C(*)*',
  };

  /** TODO: Full?
   */
  public static AAFullNames: StringDictionary = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Cys': 'C',
    'Gln': 'Q',
    'Glu': 'E',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
  };

  public static getPalette(scheme: string = 'grok'): SeqPalette {
    switch (scheme) {
    case 'grok':
      return AminoacidsPalettes.GrokGroups;
    case 'lesk':
      return AminoacidsPalettes.Lesk;
    default:
      throw new Error(`ChemPalette: scheme \`${scheme}\` does not exist`);
    }
  }

  /**
   * Returns divided amino acid with its content in the bracket, if the content is number, then its omitted
   *
   * @param {string} c raw amino
   * @return {[string, string]} outer and inner content
   */
  public static getInnerOuter(c: string): [string, string] {
    let isInner = 0;
    let inner = '';
    let outer = '';

    for (const char of c) {
      if (char == '(')
        isInner++;
      else if (char == ')')
        isInner--;
      else if (isInner)
        inner += char;
      else
        outer += char;
    }

    return !isNaN(parseInt(inner)) ? [outer, ''] : [outer, inner];
  }

  public static getColorAAPivot(monomer: string = '', scheme: 'grok' = 'grok'): [string, string, string, number] {
    //const chemPaletteInstance = AAPalettes.GrokGroups();
    const chemPaletteInstance = this.getPalette(scheme);
    let [outerMonomer, innerMonomer] = this.getInnerOuter(monomer);
    outerMonomer = (outerMonomer.length > 6 ? `${outerMonomer.slice(0, 3)}...` : outerMonomer);
    innerMonomer = (innerMonomer.length > 6 ? `${innerMonomer.slice(0, 3)}...` : innerMonomer);

    if (monomer.length == 1 || monomer[1] == '(') {
      const amino = monomer[0]?.toUpperCase()!;
      return amino in chemPaletteInstance ?
        [chemPaletteInstance.get(amino), amino, innerMonomer, 1] :
        [this.undefinedColor, outerMonomer, innerMonomer, 1];
    }

    if (monomer[0] == 'd' && monomer[1]! in chemPaletteInstance) {
      if (monomer.length == 2 || monomer[2] == '(') {
        const amino = monomer[1]?.toUpperCase()!;
        return amino in chemPaletteInstance ?
          [chemPaletteInstance.get(amino), amino, innerMonomer, 2] :
          [this.undefinedColor, outerMonomer, innerMonomer, 2];
      }
    }

    if (monomer.substring(0, 3) in this.AAFullNames) {
      if (monomer.length == 3 || monomer[3] == '(') {
        const amino = this.AAFullNames[monomer.substring(0, 3)];
        return amino in chemPaletteInstance ?
          [chemPaletteInstance.get(amino), amino, innerMonomer, 3] :
          [this.undefinedColor, outerMonomer, innerMonomer, 3];
      }
    }

    if (monomer[0]?.toLowerCase() == monomer[0]) {
      if (monomer.substring(1, 3) in this.AAFullNames) {
        if (monomer.length == 4 || monomer[4] == '(') {
          const amino = this.AAFullNames[monomer.substring(1, 3)];
          return amino in chemPaletteInstance ?
            [chemPaletteInstance.get(amino), amino, innerMonomer, 4] :
            [this.undefinedColor, outerMonomer, innerMonomer, 4];
        }
      }
    }

    return [this.undefinedColor, outerMonomer, innerMonomer, 0];
  }
}
