import * as DG from 'datagrok-api/dg';

export type PaletteCodes = {
    hex: string[],
    rgb: string[],
    rgbCut: string[],
    numerical: number[]
  };

export function getPalette(activityNum: number): PaletteCodes {
  const standradPal = DG.Color.categoricalPalette;
  const hex = Array<string>(activityNum);
  const rgb = Array<string>(activityNum);
  const rgbCut = Array<string>(activityNum);
  const numerical = Array<number>(activityNum);

  for (let i = 0; i < activityNum; i++) {
    const modNum = i % standradPal.length;
    hex[i] = DG.Color.toHtml(standradPal[modNum]);
    rgb[i] = DG.Color.toRgb(standradPal[modNum]);
    rgbCut[i] = rgb[i].replace('rgb(', '').replace(')', '');
    numerical[i] = standradPal[modNum];
  }

  return {hex, rgb, rgbCut, numerical};
}
