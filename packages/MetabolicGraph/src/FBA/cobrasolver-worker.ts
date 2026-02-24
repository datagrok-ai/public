/* eslint-disable no-undef */

import {modelFromWorkerData} from './cobra-model'
import {optimizeModel} from './cobraSolver'

onmessage = async (message) => {
  const model = modelFromWorkerData(message.data)
  try {
    console.log(message.data)
    const solution = optimizeModel(model)
    postMessage(solution)
  } catch (e) {
    postMessage({ error: e });
  }
}
