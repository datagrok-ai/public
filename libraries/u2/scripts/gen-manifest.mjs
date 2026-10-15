/* Writes manifest.json from the registry — the machine-readable component list an LLM reads
   before emitting a dg-ui/1 spec — and the kind catalog (src/core/kinds.ts) the bdd library generates
   its vocabulary from. Run `npm run build` first: this loads the compiled sources.
   With --check it compares instead of writing, so a stale checked-in file fails the build. */
import {readFileSync, realpathSync, writeFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import '../tests/dom-shim.js';

/** The catalog states a tag reaches through a prop it declares: state → [prop, the props that set it]. */
const STATES = {
  disabled: ['enabled', {enabled: false}],
  readOnly: ['readOnly', {readOnly: true}],
  invalid: ['required', {required: true}],
  busy: ['busy', {busy: true}],
  empty: ['items', {items: []}],
};

/** States a tag's fresh value keeps it out of: a checkbox holds `false`, a color input black and a
 * font input the default font, so `required` alone leaves them valid. */
const UNREACHED = {invalid: ['u2-bool-input', 'u2-color-input', 'u2-font-input']};

/** The manifest the checked-in file must equal. `base` is where the compiled sources are: `dist/src/`
 * for the script, the in-place emit (`src/`) for the drift test. */
export async function buildManifest(base = '../dist/src/') {
  const load = (path) => import(new URL(base + path, import.meta.url).href);
  const [{registry}, {registerAll}, {SpecContext, renderSpec}, {U2_KINDS, kindCatalog}] = await Promise.all([
    load('spec/registry.js'), load('spec/registrations.js'), load('spec/spec.js'), load('core/kinds.js')]);
  registerAll();
  const manifest = registry.manifest();

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

  const components = manifest.components.map((c) => {
    const names = c.props.map((p) => p.name);
    const states = Object.fromEntries(Object.entries(STATES)
      .filter(([state, [prop]]) => c.visual !== false && names.includes(prop) &&
        !(UNREACHED[state] ?? []).includes(c.tag))
      .map(([state, [, props]]) => [state, props]));
    return {...c, ...kindOfTag.has(c.tag) ? {kind: kindOfTag.get(c.tag)} : {},
      ...Object.keys(states).length > 0 ? {states} : {}};
  });
  const propsOf = (tag) => components.find((c) => c.tag === tag).props
    .filter((p) => p.structural !== true).map((p) => p.name);
  return {...manifest, components, kinds: kindCatalog(U2_KINDS, tagOf, propsOf), kindsVersion: 1};
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
    /* the checked-in copy is CRLF in a Windows checkout, the generated string always LF */
    const lf = (s) => s.replace(/\r\n/g, '\n');
    const actual = lf(readFileSync(path, 'utf8'));
    if (actual === lf(json))
      console.log(`manifest.json is up to date (${summary})`);
    else {
      const a = actual.split('\n');
      const b = lf(json).split('\n');
      let i = 0;
      while (i < a.length && i < b.length && a[i] === b[i]) i++;
      console.error('manifest.json is stale — run `npm run manifest` and commit the result.');
      console.error(`  first difference at line ${i + 1}:`);
      console.error(`    checked in: ${i < a.length ? a[i] : '<end of file>'}`);
      console.error(`    generated:  ${i < b.length ? b[i] : '<end of file>'}`);
      process.exit(1);
    }
  }
}
