import { sampleReactionsWasm } from "./sampler-wrapper";

onmessage = async (message) => {
    // split workers, start and end denote the range of reactions to be processed by this worker
    // note, the range is not inclusive and contains 0 - reactions.length * 2 to access all reactions and directions
    //
  const {model, samples, thinning} = message.data;
  try {
    const modelResult = await sampleReactionsWasm(model, samples, thinning);
    postMessage(modelResult);
  } catch (e) {
    postMessage({ error: e });
  }
}