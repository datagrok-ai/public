/* eslint-disable no-undef */

import * as cobra from './cobraSolver'

onmessage = async (message) => {
  const model = cobra.modelFromWorkerData(message.data)
  try {
    console.log(message.data)
    const solution = (model.optimize())
    postMessage(solution)
  } catch (e) {
    postMessage({ error: e });
  }
}
