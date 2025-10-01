/* eslint-disable no-undef */

import * as cobra from './glpkFBA'

onmessage = async (message) => {
  
  try {
    const modelResult = await cobra.runGLPKFBA(message.data, false);
    postMessage({fluxes: modelResult.solutions[0], reactionNames: modelResult.reactionNames});
  } catch (e) {
    postMessage({ error: e });
  }
}
