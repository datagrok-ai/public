import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SeqPalette, SeqPaletteBase} from './seq-palettes';
import {getMonomerLibHelper} from './monomer-works/monomer-utils';
import {PolymerType} from './helm/types';
import {PolymerTypes} from './helm/consts';

/** makes the color less white, makes the transparency effect always perceptible
 * @param {string} color color in string format either hex or rgb.
 * @param {boolean} scale if scale is needed to 210 brightness.
 * @return {string} color in rgb format.
 * */
function correctColor(color: string | null, scale = true): string {
  if (color == null)
    return 'rgb(100,100,100)';

  const dgColor: number = DG.Color.fromHtml(color);
  if (scale) {
    const g = DG.Color.g(dgColor);
    const r = DG.Color.r(dgColor);
    const b = DG.Color.b(dgColor);
    // calculate euclidean distance to white
    const distToBlack = Math.sqrt(Math.pow(0 - r, 2) + Math.pow(0 - g, 2) + Math.pow(0 - b, 2));
    // normalize vector r g b
    const normR = r / distToBlack;
    const normG = g / distToBlack;
    const normB = b / distToBlack;
    if (distToBlack > 210)
      return `rgb(${normR * 210},${normG * 210},${normB * 210})`;
  }

  return DG.Color.toRgb(dgColor);
}

export class StringUtils {
  public static hashCode(s: string): number {
    let hash: number = 0;
    if (s.length === 0)
      return hash;
    for (let i: number = 0; i < s.length; i++) {
      const chr: number = s.charCodeAt(i);
      hash = ((hash << 5) - hash) + chr;
      hash |= 0; // Convert to 32bit integer
    }
    return hash;
  }
}

export abstract class UnknownSeqPalette implements SeqPalette {
  public abstract get(m: string, polymerType?: string): string;
}

export class GrayAllPalette extends UnknownSeqPalette {
  public get(_m: string, _polymerType?: string): string {
    return '#666666';
  }
}


export class UnknownColorPalette extends UnknownSeqPalette {
  public static palette: string[] = UnknownColorPalette.buildPalette();
  // this way is just more future-proof, when we start distinguishing
  // between different polymer types for coloring of non natural aa's or nucleotides
  public static customMonomerColors: { [symbol: string]: { [polymerType: string]: string } } = {};
  public static polymerTypes: string[] = [];

  private static buildPalette(): string[] {
    getMonomerLibHelper().then((lh) => {
      lh.awaitLoaded(Infinity).then(() => {
        const monLib = lh.getMonomerLib();
        monLib.onChanged.subscribe(() => {
          UnknownColorPalette.customMonomerColors = {};
          UnknownColorPalette.polymerTypes = monLib.getPolymerTypes();
          for (const polymerType of this.polymerTypes) {
            const monomerSymbols = monLib.getMonomerSymbolsByType(polymerType as PolymerType);
            for (const monomerSymbol of monomerSymbols) {
              const monomer = monLib.getMonomer(polymerType as PolymerType, monomerSymbol);
              if (monomer?.meta?.colors?.default?.background) {
                if (!this.customMonomerColors[monomerSymbol])
                  this.customMonomerColors[monomerSymbol] = {};

                this.customMonomerColors[monomerSymbol][polymerType] =
                  correctColor(monomer.meta.colors.default.background);
              }
            }
          }
        });
      });
    });
    const res = ([] as string[]).concat(...Object.values(SeqPaletteBase.colourPalette));
    return res;
  }

  public get(m: string, polymerType?: string): string {
    const colorObj = UnknownColorPalette.customMonomerColors[m];
    const polType = polymerType ?? PolymerTypes.PEPTIDE;
    if (colorObj && colorObj[polType])
      return colorObj[polType];
    const hash: number = StringUtils.hashCode(m);
    const pI = hash % UnknownColorPalette.palette.length;
    return correctColor(UnknownColorPalette.palette[pI]);
  }
}

export class UnknownSeqPalettes extends SeqPaletteBase {
  private static gray: SeqPalette;

  public static get Gray(): SeqPalette {
    if (this.gray === void 0)
      this.gray = new GrayAllPalette();
    return this.gray;
  }

  private static color: SeqPalette;

  public static get Color(): SeqPalette {
    if (this.color === void 0)
      this.color = new UnknownColorPalette();
    return this.color;
  }
}
