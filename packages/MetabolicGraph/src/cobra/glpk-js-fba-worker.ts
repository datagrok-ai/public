import * as cobra from './glpkJS';

onmessage = async (message) => {
  try {
    const modelResult = await cobra.solveUsingGLPKJvail(message.data);
    postMessage({fluxes: modelResult.fluxes[0], reactionNames: modelResult.reactionNames});
  } catch (e) {
    postMessage({error: e});
  }
};
