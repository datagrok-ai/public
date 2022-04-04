import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export class ChemPalette {
    cp:{ [Key: string]: string; } = {};
    constructor(scheme: string) {
      if (scheme == 'grok')
        this.cp = ChemPalette.getDatagrok();
    }

    showTooltip(cell:DG.GridCell, x:number, y:number) {
      const s = cell.cell.value as string;
      let toDisplay = [ui.divText(s)];
      // eslint-disable-next-line no-unused-vars
      const [_c, aar, _p] = this.getColorAAPivot(s);
      if (aar in ChemPalette.AASmiles) {
        if (s in ChemPalette.AANames)
          toDisplay = [ui.divText(ChemPalette.AANames[s])];

        if (s in ChemPalette.AAFullNames)
          toDisplay = [ui.divText(ChemPalette.AANames[ChemPalette.AAFullNames[s]])];

        const sketch = grok.chem.svgMol(ChemPalette.AASmiles[aar]);
        toDisplay.push(sketch);
      }
      ui.tooltip.show(ui.divV(toDisplay), x, y);
    }

    getColor( c: string) {
    // eslint-disable-next-line no-unused-vars
      const [color, _] = this.getColorPivot(c);
      return color;
    }
    getColorAAPivot(c = ''): [string, string, number] {
      if (c.length == 1 || c.at(1) == '(') {
        const amino = c.at(0)?.toUpperCase()!;
        return amino in this.cp?
          [this.cp[amino], amino, 1]:
          [ChemPalette.undefinedColor, '', 1];
      }
      if (c.at(0) == 'd' && c.at(1)! in this.cp) {
        if (c.length == 2 || c.at(2) == '(') {
          const amino = c.at(1)?.toUpperCase()!;
          return amino in this.cp?
            [this.cp[amino], amino, 2]:
            [ChemPalette.undefinedColor, '', 2];
        }
      }
      if (c.substr(0, 3) in ChemPalette.AAFullNames) {
        if (c.length == 3 || c.at(3) == '(') {
          const amino = ChemPalette.AAFullNames[c.substr(0, 3)];
          return amino in this.cp?
            [this.cp[amino], amino, 3]:
            [ChemPalette.undefinedColor, '', 3];
        }
      }
      if (c.at(0)?.toLowerCase() == c.at(0)) {
        if (c.substr(1, 3) in ChemPalette.AAFullNames) {
          if (c.length == 4 || c.at(4) == '(') {
            const amino = ChemPalette.AAFullNames[c.substr(1, 3)];
            return amino in this.cp?
              [this.cp[amino], amino, 4]:
              [ChemPalette.undefinedColor, '', 4];
          }
        }
      }
      return [ChemPalette.undefinedColor, '', 0];
      //return c ? DG.Color.toRgb(this.colorScale(c)) : 'rgb(127,127,127)'
    }
    getColorPivot(c = ''): [string, number] {
      // eslint-disable-next-line no-unused-vars
      const [color, _, pivot] = this.getColorAAPivot(c);
      return [color, pivot];
    };

    static colourPalette: { [key: string]: string[] } = {
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
    }

    static grokGroups: [string[], string][] = [
      [['C', 'U'], 'yellow'],
      [['G', 'P'], 'red'],
      [['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'], 'all_green'],
      [['R', 'H', 'K'], 'light_blue'],
      [['D', 'E'], 'dark_blue'],
      [['S', 'T', 'N', 'Q'], 'orange'],
    ];
    static undefinedColor = 'rgb(100,100,100)'

    static makePalette(dt: [string[], string][], simplified = false) {
      const palette: { [key: string]: string } = {};
      dt.forEach((cp) => {
        const objList = cp[0];
        const colour = cp[1];
        objList.forEach((obj, ind) => {
          palette[obj] = ChemPalette.colourPalette[colour][simplified ? 0 : ind];
        });
      });
      return palette;
    }

    static AANames:{ [Key: string]: string; } = {
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
    }
    static AASmiles: { [Key: string]: string; } = {
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
      'R': 'N[C(=O)OH](CCCNC(N)=N)C(=O)O',
      'H': 'NC(CC1=CN=C[N]1)C(=O)O',
      'C': 'N[C@@H](CS)C(=O)O',
      'V': 'NC(C(C)C)C(=O)O',
      'P': 'N(CCC1)C1C(=O)O',
      'W': 'N[C@@H](Cc1c2ccccc2n([H])c1)C(=O)O',
      'I': 'N[C@H]([C@H](C)CC)C(=O)O',
      'M': 'NC(CCSC)C(=O)O',
      'T': 'NC(C(O)C)C(=O)O',
    }
    static AASmilesTruncated: { [Key: string]: string; } = {
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
      'R': '*[C(=O)OH](CCCNC(N)=N)*',
      'H': 'C1=C(NC=N1)CC(*)*',
      'C': 'C([C@@H](*)*)S',
      'V': 'CC(C)C(*)*',
      'P': 'C1CCN(*)C1*',
      'W': '*[C@@H](Cc1c2ccccc2n([H])c1)*',
      'I': 'CC[C@H](C)[C@H](*)*',
      'M': 'CSCCC(*)*',
      'T': 'CC(O)C(*)*',
    }

    static AAFullNames: { [Key: string]: string; } = {
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

    static getDatagrok = ():{ [Key: string]: string; } => {
      return ChemPalette.makePalette(ChemPalette.grokGroups);
    };

    getLesk = () => ChemPalette.makePalette([
      [['G', 'A', 'S', 'T'], 'orange'],
      [['C', 'V', 'I', 'L', 'P', 'F', 'Y', 'M', 'W'], 'all_green'],
      [['N', 'Q', 'H'], 'magenta'],
      [['D', 'E'], 'red'],
      [['K', 'R'], 'all_blue']])
}
