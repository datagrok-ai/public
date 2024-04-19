import {optimizeNM} from './optimizer-nelder-mead';
import {Extremum} from './optimizer-misc';

onmessage = async (event) => {
  const {objectiveFunc, params, start, end} = event.data;
  const extremums = new Array<Extremum>(end - start);
  const warnings = new Array<string>(0);
  for (let parametrization = start; parametrization < end; ++parametrization) {
    try {
      const extremum = await optimizeNM(objectiveFunc, params[parametrization], {});
      extremums[parametrization - start] = extremum;
    } catch (err: any) {
      const errMsg: string = err instanceof Error ? err.message : err.toString();
      warnings.push(errMsg);
    }
  }
  postMessage({extremums, warnings});
};
