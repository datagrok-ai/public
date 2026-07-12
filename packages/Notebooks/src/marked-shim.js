// marked v4 (security override for GHSA-rrrm-qjm4-v8hf / GHSA-5v2h-r2cx-5xgj) removed the
// default export that @jupyterlab/rendermime@2.x relies on (`import marked from 'marked'`).
// Re-export the callable `marked` function as default; it still carries setOptions/Renderer/etc.
// Wired up via `resolve.alias` in webpack.config.js.
import {marked} from 'marked-esm';

export default marked;
export * from 'marked-esm';
