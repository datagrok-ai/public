// ESM module-customization hook used by the Node test runner (package-test-node.ts).
// The js-api / test sources are written for the browser and statically import CSS
// (and occasionally other assets), which webpack handles via loaders but Node cannot
// parse. This hook returns an empty module for any asset import so the runner can load
// the same TypeScript sources the browser bundles, unchanged.
const ASSET = /\.(css|scss|sass|less|svg|png|jpe?g|gif|webp|woff2?|ttf|eot)(\?.*)?$/;

export async function load(url, context, nextLoad) {
  if (ASSET.test(url))
    return {format: 'module', source: 'export default {};', shortCircuit: true};
  return nextLoad(url, context);
}
