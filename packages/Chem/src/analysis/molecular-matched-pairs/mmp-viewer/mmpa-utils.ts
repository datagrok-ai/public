import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SCALING_METHODS} from './mmp-constants';

/**
 * Scales activity column values.
 * @param activityCol - Activity column.
 * @param scaling - Scaling method.
 * @return - Scaled activity column.
 */
export function scaleActivity(activityCol: DG.Column<number>, scaling: SCALING_METHODS = SCALING_METHODS.NONE,
): DG.Column<number> {
  let formula = (x: number): number => x;
  switch (scaling) {
  case SCALING_METHODS.NONE:
    break;
  case SCALING_METHODS.LG:
    formula = (x: number): number => Math.log10(x);
    break;
  case SCALING_METHODS.MINUS_LG:
    formula = (x: number): number => -Math.log10(x);
    break;
  default:
    throw new Error(`ScalingError: method \`${scaling}\` is not available.`);
  }
  const activityColData = activityCol.getRawData();
  const scaledCol: DG.Column<number> = DG.Column.float(`${scaling} ${activityCol.name}`, activityCol.length)
    .init((i) => {
      const val = activityColData[i];
      return val === DG.FLOAT_NULL || val === DG.INT_NULL ? val : formula(val);
    });
  scaledCol.setTag(DG.TAGS.FORMULA, scaling);
  return scaledCol;
}
