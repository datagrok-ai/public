import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {MonomerLibrary} from '../monomer-library';


export class ChemPalette {
  cp: StringDictionary = {};
  isInit: boolean = false;
  monomerLib: MonomerLibrary | null = null;

  /**
   * Creates an instance of ChemPalette.
   *
   * @param {string} scheme Color scheme to use.
   * @param {boolean} [grouping=false] Is grouping enabled.
   * @memberof ChemPalette
   */
  private constructor(scheme: string, grouping = false) {
    if (scheme == 'grok')
      this.cp = ChemPalette.getDatagrok(grouping);
  }

  /**
   * Renders 2D representation of a amino acid residue in a tooltip.
   *
   * @param {DG.GridCell} cell Grid cell to show tooltip over.
   * @param {number} x x coordinate of the mouse pointer.
   * @param {number} y y coordinate of the mouse pointer.
   */
  static showTooltip(cell: DG.GridCell, x: number, y: number, monomerLib: MonomerLibrary) {
    // if (!this.isInit) {
    //   const validPackage = _package ?? _packageTest;
    //   if (!validPackage)
    //     throw new Error('No package instance found');
    //   this.monomerLib = new MonomerLibrary(await validPackage.files.readAsText(`HELMMonomers_June10.sdf`));
    //   this.isInit = true;
    // }

    const s = cell.cell.value as string;
    let toDisplay = [ui.divText(s)];
    const [, aarOuter, aarInner] = ChemPalette.getColorAAPivot(s);
    for (const aar of [aarOuter, aarInner]) {
      if (monomerLib.monomerNames.includes(aar)) {
        if (aar in ChemPalette.AANames)
          toDisplay = [ui.divText(ChemPalette.AANames[aar])];

        if (aar in ChemPalette.AAFullNames)
          toDisplay = [ui.divText(ChemPalette.AANames[ChemPalette.AAFullNames[aar]])];

        const options = {
          autoCrop: true,
          autoCropMargin: 0,
          suppressChiralText: true,
        };
        const sketch = grok.chem.svgMol(monomerLib.getMonomerMol(aar), undefined, undefined, options);
        if (toDisplay.length == 2)
          toDisplay.push(ui.divText('Modified'));

        toDisplay.push(sketch);
      }
    }
    ui.tooltip.show(ui.divV(toDisplay), x, y);
  }

  /**
   * Retursn divided amino with its content in the bracket, if the conetent is number, then its omitted
   *
   * @param {string} c raw amino
   * @return {[string, string]} outer and inner content
   */
  static getInnerOuter(c: string): [string, string] {
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

  static getColorAAPivot(monomer: string = '', scheme: 'grok' = 'grok'): [string, string, string, number] {
    // const chemPaletteInstance = ChemPalette.getDatagrok();
    const chemPaletteInstance = ChemPalette.getPalette(scheme);
    let [outerMonomer, innerMonomer] = ChemPalette.getInnerOuter(monomer);
    outerMonomer = (outerMonomer.length > 6 ? `${outerMonomer.slice(0, 3)}...` : outerMonomer);
    innerMonomer = (innerMonomer.length > 6 ? `${innerMonomer.slice(0, 3)}...` : innerMonomer);

    if (monomer.length == 1 || monomer[1] == '(') {
      const amino = monomer[0]?.toUpperCase()!;
      return amino in chemPaletteInstance ?
        [chemPaletteInstance[amino], amino, innerMonomer, 1]:
        [ChemPalette.undefinedColor, outerMonomer, innerMonomer, 1];
    }

    if (monomer[0] == 'd' && monomer[1]! in chemPaletteInstance) {
      if (monomer.length == 2 || monomer[2] == '(') {
        const amino = monomer[1]?.toUpperCase()!;
        return amino in chemPaletteInstance ?
          [chemPaletteInstance[amino], amino, innerMonomer, 2]:
          [ChemPalette.undefinedColor, outerMonomer, innerMonomer, 2];
      }
    }

    if (monomer.substring(0, 3) in ChemPalette.AAFullNames) {
      if (monomer.length == 3 || monomer[3] == '(') {
        const amino = ChemPalette.AAFullNames[monomer.substring(0, 3)];
        return amino in chemPaletteInstance ?
          [chemPaletteInstance[amino], amino, innerMonomer, 3]:
          [ChemPalette.undefinedColor, outerMonomer, innerMonomer, 3];
      }
    }

    if (monomer[0]?.toLowerCase() == monomer[0]) {
      if (monomer.substring(1, 3) in ChemPalette.AAFullNames) {
        if (monomer.length == 4 || monomer[4] == '(') {
          const amino = ChemPalette.AAFullNames[monomer.substring(1, 3)];
          return amino in chemPaletteInstance ?
            [chemPaletteInstance[amino], amino, innerMonomer, 4]:
            [ChemPalette.undefinedColor, outerMonomer, innerMonomer, 4];
        }
      }
    }

    return [ChemPalette.undefinedColor, outerMonomer, innerMonomer, 0];
  }

  static colourPalette: {[key: string]: string[]} = {
    'orange': ['rgb(255,187,120)', 'rgb(245,167,100)', 'rgb(235,137,70)', 'rgb(205, 111, 71)'],
    'all_green': ['rgb(44,160,44)', 'rgb(74,160,74)', 'rgb(23,103,57)', 'rgb(30,110,96)', 'rgb(60,131,95)',
      'rgb(24,110,79)', 'rgb(152,223,138)', 'rgb(182, 223, 138)', 'rgb(152, 193, 138)'],
    'all_blue': ['rgb(31,119,180)', 'rgb(23,190,207)', 'rgb(122, 102, 189)', 'rgb(158,218,229)', 'rgb(141, 124, 217)',
      'rgb(31, 120, 150)'],
    'magenta': ['rgb(162,106,192)', 'rgb(197,165,224)', 'rgb(208,113,218)'],
    'red': ['rgb(214,39,40)', 'rgb(255,152,150)'],
    'st_blue': ['rgb(23,190,207)', 'rgb(158,218,229)', 'rgb(31,119,180)'],
    'dark_blue': ['rgb(31,119,180)', 'rgb(31, 120, 150)'],
    'light_blue': ['rgb(23,190,207)', 'rgb(158,218,229)', 'rgb(108, 218, 229)', 'rgb(23,190,227)'],
    'lilac_blue': ['rgb(124,102,211)', 'rgb(149,134,217)', 'rgb(97, 81, 150)'],
    'dark_green': ['rgb(23,103,57)', 'rgb(30,110,96)', 'rgb(60,131,95)', 'rgb(24,110,79)'],
    'green': ['rgb(44,160,44)', 'rgb(74,160,74)'],
    'light_green': ['rgb(152,223,138)', 'rgb(182, 223, 138)', 'rgb(152, 193, 138)'],
    'st_green': ['rgb(44,160,44)', 'rgb(152,223,138)', 'rgb(39, 174, 96)', 'rgb(74,160,74)'],
    'pink': ['rgb(247,182,210)'],
    'brown': ['rgb(140,86,75)', 'rgb(102, 62, 54)'],
    'gray': ['rgb(127,127,127)', 'rgb(199,199,199)', 'rgb(196,156,148)', 'rgb(222, 222, 180)'],
    'yellow': ['rgb(188,189,34)'],
    'white': ['rgb(230,230,230)'],
  };

  static grokGroups: {[key: string]: string[]} = {
    'yellow': ['C', 'U'],
    'red': ['G', 'P'],
    'all_green': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
    'light_blue': ['R', 'H', 'K'],
    'dark_blue': ['D', 'E'],
    'orange': ['S', 'T', 'N', 'Q'],
  };

  static undefinedColor = 'rgb(100,100,100)';

  static makePalette(dt: {[key: string]: string[]}, simplified = false, grouping = false): StringDictionary {
    const palette: { [key: string]: string } = {};
    const groups = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    let currentGroup = 0;
    for (const [color, monomers] of Object.entries(dt)) {
      monomers.forEach((monomer, index) => {
        palette[grouping ? groups[currentGroup] : monomer] = ChemPalette.colourPalette[color][simplified ? 0 : index];
      });
      currentGroup++;
    }
    return palette;
  }

  static AANames: StringDictionary = {
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

  static AASmiles: StringDictionary = {
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

  static AASmilesTruncated: StringDictionary = {
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

  static AAFullNames: StringDictionary = {
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

  static getDatagrok(grouping = false): StringDictionary {
    return ChemPalette.makePalette(ChemPalette.grokGroups, false, grouping);
  }

  static getPalette(scheme: 'grok'): StringDictionary {
    switch (scheme) {
    case 'grok':
      return ChemPalette.getDatagrok();
    default:
      throw new Error(`ChemPalette: scheme \`${scheme}\` does not exist`);
    }
  }
}
