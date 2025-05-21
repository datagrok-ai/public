// Worker for sensitivity analysis of Diff Studio model

import {NO_ERRORS, RESULT_CODE, WorkerTask} from '../defs';
import {getAnalysis} from '../sens-analysis-utils';


onmessage = async function(evt) {
  try {
    const task = evt.data.task as WorkerTask;
    const batchSize = task.pipelines.length;
    const rowVals = new Array<Float64Array>(batchSize);
    const analysisRes = new Array<string>(batchSize);

    for (let idx = 0; idx < batchSize; ++idx) {
      try {
        rowVals[idx] = getAnalysis({
          ivp2ww: task.ivp2ww,
          pipeline: task.pipelines[idx],
          inputVec: task.inputVecs[idx],
          rowIdx: task.rowIdx,
          argColIdx: task.argColIdx,
          argVal: task.argVal,
        });

        analysisRes[idx] = NO_ERRORS;
      } catch (e) {
        analysisRes[idx] = e instanceof Error ? e.message : 'Platform issue: in-webworker analysis failed';
      }
    }

    postMessage({
      'callResult': RESULT_CODE.SUCCEED,
      'rowVals': rowVals,
      'analysisRes': analysisRes,
    });
  } catch (e) {
    postMessage({
      'callResult': RESULT_CODE.FAILED,
      'msg': e instanceof Error ? e.message : 'Platform issue: in-webworker analysis failed',
    });
  }
};
