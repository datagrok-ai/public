import {Extremum} from '../../optimizer-misc';
import {fit} from '../fitting-utils';
import {NO_ERRORS, NelderMeadInput} from '../defs';

onmessage = async function(evt) {
  try {
    const task = evt.data.task as NelderMeadInput;
    const startPoints = evt.data.startPoints as Float32Array[];
    const pointsCount = startPoints.length;
    const fitRes = new Array<string>(pointsCount);
    const extremums: Extremum[] = [];

    for (let i = 0; i < pointsCount; ++i) {
      try {
        extremums.push(await fit(task, startPoints[i]));
        fitRes[i] = NO_ERRORS;
      } catch (e) {
        fitRes[i] = e instanceof Error ? e.message : 'Platform issue: in-webworker optimization failed';
      }
    }

    postMessage({'callResult': 0, 'extremums': extremums, 'fitRes': fitRes});
  } catch (e) {
    postMessage({
      'callResult': -1,
      'msg': e instanceof Error ? e.message : 'Platform issue: in-webworker optimization failed',
    });
  }
};
