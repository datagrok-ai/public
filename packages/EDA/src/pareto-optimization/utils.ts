import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {OPT_TYPE} from './defs';

export const PALETTE = [DG.Color.darkGreen, DG.Color.yellow, DG.Color.darkRed];

/** Return output color palette w.r.t. the specified type of optimization */
export function getOutputPalette(type: OPT_TYPE): number[] {
  if (type === OPT_TYPE.MIN)
    return [...PALETTE];

  return [...PALETTE].reverse();
}

/** Return div with color scale description */
export function getColorScaleDiv(type: OPT_TYPE, useMinMax: boolean = true): HTMLElement {
  const scale = ui.label('Color scale:');
  scale.style.paddingRight = '7px';
  const elems = [scale];
  const minLbl = ui.label(useMinMax ? 'min' : 'worst');
  const midLbl = ui.label('. . .');
  const maxLbl = ui.label(useMinMax ? 'max' : 'best');
  const palette = getOutputPalette(type);

  const colorElems = [minLbl, midLbl, maxLbl].map((el, idx) => {
    if (idx !== 1) {
      el.style.fontWeight = 'bold';
      el.style.color = DG.Color.toRgb(palette[idx]);
    }

    el.style.marginRight = '5px';

    return el;
  });

  elems.push(...colorElems);

  return ui.divH(elems);
} // getColorScaleDiv
