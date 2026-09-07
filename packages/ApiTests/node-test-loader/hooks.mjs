// ESM module-customization hook used by the Node test runner (package-test-node.ts).
// The js-api / test sources are written for the browser and statically import CSS
// (and occasionally other assets), which webpack handles via loaders but Node cannot
// parse. This hook returns an empty module for any asset import so the runner can load
// the same TypeScript sources the browser bundles, unchanged.
//
// The datagrok-api/{dg,grok,ui} aliasing is NOT here: hooks run on their own thread and
// cannot read the runtime globals, and whether this hook or tsx's resolves a given import
// is not something to depend on. register.mjs binds those to generated modules instead.
const ASSET = /\.(css|scss|sass|less|svg|png|jpe?g|gif|webp|woff2?|ttf|eot)(\?.*)?$/;

export async function load(url, context, nextLoad) {
  if (ASSET.test(url))
    return {format: 'module', source: 'export default {};', shortCircuit: true};
  return nextLoad(url, context);
}
