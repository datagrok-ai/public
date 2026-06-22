// Stubs browser asset imports so the Node test runner can load the same browser-oriented
// js-api / test TypeScript sources webpack bundles, unchanged. Loaded via
// `node --import tsx --import ./node-test-loader/register.mjs` (after tsx, so the CJS
// extension overrides below win over tsx's transformer).
import {register} from 'node:module';
import Module from 'node:module';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc.js';
import advancedFormat from 'dayjs/plugin/advancedFormat.js';

// ESM side: intercept `import './x.css'` from ESM-graph modules.
register('./hooks.mjs', import.meta.url);

// Browser env shim: the test/js-api sources assume dayjs plugins are already loaded.
dayjs.extend(utc);
dayjs.extend(advancedFormat);

// CJS side: intercept `require('./x.css')` from modules tsx treats as CommonJS.
const ASSET_EXTENSIONS = [
  '.css', '.scss', '.sass', '.less',
  '.svg', '.png', '.jpg', '.jpeg', '.gif', '.webp',
  '.woff', '.woff2', '.ttf', '.eot',
];
const stub = (module) => { module.exports = {}; };
for (const ext of ASSET_EXTENSIONS)
  Module._extensions[ext] = stub;
