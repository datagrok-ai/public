import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

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
    return UnknownColorPalette.palette[pI];
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
