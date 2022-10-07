import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** makes the color less white, makes the transparency effect always perceptible
 * @param {string} color x coordinate.
 * */
function correctColor(color: string | null): string {
  if (color == null)
    return 'rgb(100,100,100)';

  const dgColor: number = DG.Color.fromHtml(color);
  const g = DG.Color.g(dgColor);
  const r = DG.Color.r(dgColor);
  const b = DG.Color.b(dgColor);
  // calculate euclidean distance to white
  const distToBlack = Math.sqrt(Math.pow(0 - r, 2) + Math.pow(0 - g, 2) + Math.pow(0 - b, 2));
  // normalize vector r g b
  const normR = r / distToBlack;
  const normG = g / distToBlack;
  const normB = b / distToBlack;
  if (distToBlack > 210) {
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

import {SeqPalette, SeqPaletteBase} from './seq-palettes';

export abstract class UnknownSeqPalette implements SeqPalette {
  public abstract get(m: string): string;
}

export class GrayAllPalette extends UnknownSeqPalette {
  public get(m: string): string {
    return '#666666';
  }
}


export class UnknownColorPalette extends UnknownSeqPalette {
  public static palette: string[] = UnknownColorPalette.buildPalette();

  private static buildPalette(): string[] {
    const res = ([] as string[]).concat(...Object.values(SeqPaletteBase.colourPalette));
    return res;
  }

  public get(m: string): string {
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
