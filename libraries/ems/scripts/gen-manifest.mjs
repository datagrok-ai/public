// EMS's components and kinds only; the format's preamble is in u2's manifest. `--check` fails on a stale file.
import {readFileSync, realpathSync, writeFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import '@datagrok-libraries/u2/tests/dom-shim.js';
import {Registry, SpecContext, renderSpec, kindCatalog} from '@datagrok-libraries/u2';

/** The manifest the checked-in file must equal. `base` is where the compiled sources are: `dist/src/`
 * for the script, the in-place emit (`src/`) for the drift test. */
export async function buildManifest(base = '../dist/src/') {
  const load = (path) => import(new URL(base + path, import.meta.url).href);
  const [{registerDomainSources}, {EMS_KINDS}] = await Promise.all([
    load('core/registrations.js'), load('core/kinds.js')]);
  const registry = new Registry();
  registerDomainSources(registry);
  const {components} = registry.manifest();

  // a tag's kind is the data-u2 its built example root carries
  const kindOfTag = new Map();
  for (const meta of registry.metas()) {
    if (meta.visual === false)
      continue;
    const instance = renderSpec({$schema: 'dg-ui/1', root: meta.example}, new SpecContext(), registry);
    const kind = instance.root.children[0]?.dataset.u2;
    if (kind !== undefined)
      kindOfTag.set(meta.tag, kind);
    instance.dispose();
  }
  const tagOf = (kind) => [...kindOfTag].find(([, k]) => k === kind)?.[0];
  const propsOf = (tag) => components.find((c) => c.tag === tag).props
    .filter((p) => p.structural !== true).map((p) => p.name);
  return {
    components: components.map((c) => kindOfTag.has(c.tag) ? {...c, kind: kindOfTag.get(c.tag)} : c),
    kinds: kindCatalog(EMS_KINDS, tagOf, propsOf),
    kindsVersion: 1,
  };
}

const isMain = process.argv[1] !== undefined &&
  realpathSync(process.argv[1]) === realpathSync(fileURLToPath(import.meta.url));

if (isMain) {
  const manifest = await buildManifest();
  const path = fileURLToPath(new URL('../manifest.json', import.meta.url));
  const json = JSON.stringify(manifest, null, 2) + '\n';
  const summary = `${manifest.components.length} components, ${manifest.kinds.length} kinds`;

  if (!process.argv.includes('--check')) {
    writeFileSync(path, json);
    console.log(`${summary} → ${path}`);
  } else {
    const lf = (s) => s.replace(/\r\n/g, '\n');
    const actual = lf(readFileSync(path, 'utf8'));
    if (actual === lf(json))
      console.log(`manifest.json is up to date (${summary})`);
    else {
      console.error('manifest.json is stale — run `npm run manifest` and commit the result.');
      process.exit(1);
    }
  }
}
