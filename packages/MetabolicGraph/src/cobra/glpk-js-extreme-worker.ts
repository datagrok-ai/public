import * as cobra from './glpkJS';

onmessage = async (message) => {
  try {
    // split workers, start and end denote the range of reactions to be processed by this worker
  // note, the range is not inclusive and contains 0 - reactions.length * 2 to access all reactions and directions
  //
    const {model, start, end} = message.data;
    const modelResult = await cobra.extremeSolveUsingGLPKJvail(model, start, end);
    postMessage(modelResult);
  } catch (e) {
    postMessage({error: e});
  }
};
