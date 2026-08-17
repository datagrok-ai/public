/* Writes manifest.json from the registry — the machine-readable component list an LLM reads
   before emitting a dg-ui/1 spec. Run `npm run build` first: this loads the compiled sources. */
import {writeFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import '../tests/dom-shim.js';
import {registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';

registerAll();
const manifest = registry.manifest();
const path = fileURLToPath(new URL('../manifest.json', import.meta.url));
writeFileSync(path, JSON.stringify(manifest, null, 2) + '\n');
console.log(`${manifest.components.length} components → ${path}`);
