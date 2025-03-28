import { FeatureExtractionPipeline, pipeline } from "@xenova/transformers";

let extractor: FeatureExtractionPipeline;

async function getExtractor() {
  if (!extractor) {
    extractor = await pipeline("feature-extraction", "Xenova/all-MiniLM-L6-v2");
  }
  return extractor;
}

async function extractFeatures(sentences: string[]) {
  const extractor = await getExtractor();
  return (await extractor(sentences, { pooling: "mean", normalize: true })).tolist();
}

self.onmessage = async ({ data: { sentences } }) => {
  let response: { error?: any; embedding?: any };
  try {
    const embedding = await extractFeatures(sentences);
    response = { embedding };
  } catch (e: any) {
    response = { error: e.toString() };
  }

  self.postMessage(response);
};