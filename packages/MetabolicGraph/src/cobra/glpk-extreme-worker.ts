/* eslint-disable no-undef */

import * as cobra from './glpkFBA'

onmessage = async (message) => {
    // split workers, start and end denote the range of reactions to be processed by this worker
    // note, the range is not inclusive and contains 0 - reactions.length * 2 to access all reactions and directions
    //
  const {model, start, end} = message.data;
  try {
    const modelResult = await cobra.runGLPKFBA(model, true, start, end);
    postMessage({fluxes: modelResult.solutions, reactionNames: modelResult.reactionNames});
  } catch (e) {
    postMessage({ error: e });
  }
}
