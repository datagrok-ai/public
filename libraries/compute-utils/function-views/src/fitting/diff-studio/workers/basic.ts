// Worker for optimization of parameters of Diff Studio model

import {Extremum} from '../../optimizer-misc';
import {fit} from '../fitting-utils';
import {NO_ERRORS, NelderMeadInput, RESULT_CODE} from '../defs';
import {COST_FUNC_THRESH, STOP_AFTER_DEFAULT} from '../../constants';

onmessage = async function(evt) {
  try {
    const task = evt.data.task as NelderMeadInput;
    const startPoints = evt.data.startPoints as Float64Array[];
    const pointsCount = startPoints.length;
    const fitRes = new Array<string>(pointsCount);
    const extremums: Extremum[] = [];
    let validPoints = 0;
    const checkValidPoints = task.earlyStoppingSettings.useEarlyStopping;
    const maxValidPoints = task.earlyStoppingSettings.stopAfter ?? STOP_AFTER_DEFAULT;
    const threshold = task.earlyStoppingSettings.costFuncThreshold ?? COST_FUNC_THRESH;

    for (let i = 0; i < pointsCount; ++i) {
      try {
        const extremum = await fit(task, startPoints[i]);
        extremums.push(extremum);
        fitRes[i] = NO_ERRORS;

        if (checkValidPoints) {
          if (extremum.cost <= threshold) {
            ++validPoints;

            if (validPoints >= maxValidPoints)
              break;
          }
        }
      } catch (e) {
        fitRes[i] = e instanceof Error ? e.message : 'Platform issue: in-webworker optimization failed';
      }
    }

    postMessage({'callResult': RESULT_CODE.SUCCEED, 'extremums': extremums, 'fitRes': fitRes});
  } catch (e) {
    postMessage({
      'callResult': RESULT_CODE.FAILED,
      'msg': e instanceof Error ? e.message : 'Platform issue: in-webworker optimization failed',
    });
  }
};
