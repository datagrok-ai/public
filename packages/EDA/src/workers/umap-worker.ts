// Worker for the method UMAP

import {UMAP} from 'umap-js';

onmessage = async function(evt) {
  const umap = new UMAP(evt.data.options);
  const embeddings = umap.fit(evt.data.data);
  postMessage({'embeddings': embeddings});
};
