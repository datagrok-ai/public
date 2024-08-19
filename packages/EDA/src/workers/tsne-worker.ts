// Worker for the method t-SNE

import {TSNE} from '@keckelt/tsne';

onmessage = async function(evt) {
  const tsne = new TSNE({
    epsilon: evt.data.options.learningRate,
    perplexity: evt.data.options.perplexity,
    dim: evt.data.options.components,
  });

  tsne.initDataRaw(evt.data.data);

  const iterCount = evt.data.options.iterations;

  for (let i = 0; i < iterCount; ++i)
    tsne.step();

  postMessage({'embeddings': tsne.getSolution()});
};
