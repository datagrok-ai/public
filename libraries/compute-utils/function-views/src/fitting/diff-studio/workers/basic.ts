import {fit, NelderMeadInput} from '../nelder-mead';

onmessage = async function(evt) {
  try {
    const task = evt.data.task as NelderMeadInput;
    const extremum = await fit(task);
    postMessage({'callResult': 0, 'res': extremum});
  } catch (e) {
    postMessage({'callResult': -1, 'msg': e instanceof Error ? e.message : ':((('});
  }
};
