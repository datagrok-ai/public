import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface SeqPalette {
  // There are too much problem with indexer property in typescript.
  // /**
  //  * @param {string} m Monomer character
  //  * @return {string} Color
  //  */
  // [m: string]: string;
  /** Monomer color
   * @param {string} m Monomer
   */
  get(m: string, polymerType?: string): string;
}

export class SeqPaletteBase implements SeqPalette {
  public static undefinedColor = 'rgb(100,100,100)';

  /** Palette with shades of primary colors */
  public static colourPalette: { [key: string]: string[] } = {
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

  protected static makePalette(dt: [string[], string][],
    simplified = false, PaletteType: typeof SeqPaletteBase = SeqPaletteBase
  ) {
    const palette: { [key: string]: string } = {};
    dt.forEach((cp) => {
      const objList = cp[0];
      const colour = cp[1];
      objList.forEach((obj, ind) => {
        palette[obj] = this.colourPalette[colour][simplified ? 0 : ind];
      });
    });
    return new PaletteType(palette);
  }

  private readonly _palette: { [m: string]: string };

  constructor(palette: { [m: string]: string }) {
    this._palette = palette;
  }

  public get(m: string, polymerType?: string): string {
    return this._palette[m];
  }
}
