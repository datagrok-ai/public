/// `grok kg check` / `grok kg gen` against the mini monorepo under fixtures/kg/good and the
/// deliberately broken home documents under fixtures/kg/broken (CONVENTIONS.md §5, §7, §10).
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {parseMember, loadTypeSystem, TypeSystem} from '../utils/kg/types';
import {splitFrontmatter} from '../utils/kg/frontmatter';
import {loadHomes, makeReport, HomeSet} from '../utils/kg/homes';
import {generate, generateDts, spliceGlossary, generateFeatures, writeOutputs} from '../utils/kg/gen';
import {kg} from '../commands/kg';

const fixtures = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures', 'kg');
const KG_DIR = path.join('core', 'docs', 'knowledge-graph');
const TYPES = ['feature', 'ticket', 'person', 'team', 'customer'];

/** A throwaway copy of the good monorepo, so a test can break one file without touching the fixture. */
function makeRepo(): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-kg-'));
  fs.cpSync(path.join(fixtures, 'good'), dir, {recursive: true});
  return dir;
}

function load(repo: string): {system: TypeSystem, homes: HomeSet} {
  const system = loadTypeSystem(path.join(repo, KG_DIR));
  return {system, homes: loadHomes(system, repo)};
}

/** Error messages of the good repo plus one broken home document dropped in at core/docs/broken.md. */
function brokenErrors(name: string): string[] {
  const repo = makeRepo();
  fs.copyFileSync(path.join(fixtures, 'broken', name), path.join(repo, 'core', 'docs', 'broken.md'));
  const {system, homes} = load(repo);
  return [...system.errors, ...homes.errors].map((e) => `${e.file}:${e.line}: ${e.message}`);
}

/** Type-system error messages after [mutate] rewrote type files in a copy of the good repo. */
function typeErrors(mutate: (kgRoot: string) => void): string[] {
  const repo = makeRepo();
  mutate(path.join(repo, KG_DIR));
  return loadTypeSystem(path.join(repo, KG_DIR)).errors.map((e) => e.message);
}

function append(kgRoot: string, file: string, text: string): void {
  fs.appendFileSync(path.join(kgRoot, file), text);
}

function write(root: string, file: string, text: string): void {
  fs.mkdirSync(path.dirname(path.join(root, file)), {recursive: true});
  fs.writeFileSync(path.join(root, file), text);
}

/** Runs the command with console captured; resets the exit code it may have set. */
async function run(argv: Record<string, unknown>): Promise<{ok: boolean, out: string[], err: string[], exitCode: number | undefined}> {
  const log = vi.spyOn(console, 'log').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  const before = process.exitCode;
  try {
    const ok = await kg(argv);
    return {ok, out: log.mock.calls.map((c) => String(c[0])), err: error.mock.calls.map((c) => String(c[0])), exitCode: process.exitCode as number | undefined};
  } finally {
    process.exitCode = before;
    log.mockRestore();
    error.mockRestore();
  }
}

describe('kg member syntax (CONVENTIONS §7.1)', () => {
  it('reads a nullable scalar', () => {
    const {member} = parseMember('slack?', 'string');
    expect(member).toMatchObject({name: 'slack', nullable: true, list: false, kind: 'scalar', scalar: 'string'});
  });

  it('reads a default on a nullable member', () => {
    expect(parseMember('count?', 'number = 0').member).toMatchObject({default: 0, scalar: 'number'});
    expect(parseMember('admin?', 'boolean = false').member).toMatchObject({default: false});
  });

  it('reads an enum with a default that is one of its literals', () => {
    const {member} = parseMember('tier?', "'standard' | 'trial' = 'standard'");
    expect(member).toMatchObject({kind: 'enum', literals: ['standard', 'trial'], default: 'standard'});
  });

  it('reads a reference union as kebab type names', () => {
    const {member} = parseMember('owner?', 'Person | Team', {types: TYPES});
    expect(member).toMatchObject({kind: 'ref', refs: ['person', 'team']});
  });

  it('reads a list', () => {
    expect(parseMember('areas?', 'Feature[]', {types: TYPES}).member).toMatchObject({kind: 'ref', refs: ['feature'], list: true});
    expect(parseMember('tags?', 'string[]').member).toMatchObject({kind: 'scalar', scalar: 'string', list: true});
  });

  it('reads exactly Record<string, string>', () => {
    expect(parseMember('meta?', 'Record<string, string>').member).toMatchObject({kind: 'record', recordValue: 'string'});
    expect(parseMember('m?', 'Record<string, Member>').error).toMatch(/only Record<string, string>/);
    expect(parseMember('m?', 'Record<string, Member>', {schema: true}).member).toMatchObject({recordValue: 'Member'});
  });

  it('rejects everything outside the vocabulary', () => {
    expect(parseMember('f', '() => void').error).toMatch(/not allowed/);
    expect(parseMember('t', '[string, number]').error).toMatch(/not allowed/);
    expect(parseMember('i', 'Person & Team', {types: TYPES}).error).toMatch(/not allowed/);
    expect(parseMember('g', 'Map<string, string>').error).toMatch(/only Record<string, string>|generic/);
    expect(parseMember('a', 'any').error).toMatch(/not allowed/);
    expect(parseMember('n', 'string[][]').error).toMatch(/nested arrays/);
    expect(parseMember('u', 'Unknown').error).toMatch(/unknown type 'Unknown'/);
    expect(parseMember('m', 'string | number').error).toMatch(/all string literals or all node types/);
  });

  it('rejects a default on a required member, a non-literal default, and a default outside the enum', () => {
    expect(parseMember('s', "string = 'x'").error).toMatch(/only on a nullable member/);
    expect(parseMember('s?', 'string = foo()').error).toMatch(/string, number or boolean literal/);
    expect(parseMember('t?', "'a' | 'b' = 'c'").error).toMatch(/'c' is not one of 'a' \| 'b'/);
    expect(parseMember('n?', "number = 'x'").error).toMatch(/is not a number/);
  });

  it('rejects a bad member name and a non-string spec', () => {
    expect(parseMember('Bad-Name', 'string').error).toMatch(/bad member name/);
    expect(parseMember('n', 5).error).toMatch(/must be a string/);
  });
});

describe('kg type files (CONVENTIONS §7.3, §7.4, schema.yaml constraints)', () => {
  it('validates the fixture type system clean and merges members down the chain', () => {
    const system = loadTypeSystem(path.join(fixtures, 'good', KG_DIR));
    expect(system.errors).toEqual([]);
    const developer = system.nodes.get('developer')!;
    expect(developer.chain).toEqual(['developer', 'person', 'actor', 'node']);
    expect(developer.prefix).toBe('P');
    expect(developer.authored).toBe(true);
    expect(developer.members.company).toMatchObject({nullable: false, refs: ['team']});
    expect(developer.members.email).toMatchObject({nullable: false});
    expect(system.edges.get('COVERS')!.properties).toHaveProperty('strength');
    expect(system.keys.get('superseded_by')!.keySide).toBe('to');
  });

  it('accepts narrowing: nullable to required, enum subset, reference to a subtype', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/developer.yaml', "  tier:          \"'staff'\"\n"))).toEqual([]);
  });

  it('rejects a widened enum', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/developer.yaml', "  tier?:         \"'staff' | 'alien'\"\n")))
      .toEqual([expect.stringMatching(/developer\.tier: adds 'alien' not allowed by the ancestor/)]);
  });

  it('rejects loosening required to nullable', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/developer.yaml', '  email?:        string\n')))
      .toEqual([expect.stringMatching(/developer\.email: loosens a required member to nullable/)]);
  });

  it('rejects a scalar kind change', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/developer.yaml', '  slack?:        number\n')))
      .toEqual([expect.stringMatching(/developer\.slack: changes string to number/)]);
  });

  it('checks a redeclared member against every ancestor, not only the parent', () => {
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'nodes/lead-developer.yaml'), [
      'type: lead-developer', 'extends: developer', 'description: Leads a team.', 'properties:',
      '  company:       Customer', '  slack?:        number', '',
    ].join('\n')))).toEqual([
      expect.stringMatching(/lead-developer\.company: references Customer, not a subtype of Team \(declared by developer/),
      expect.stringMatching(/lead-developer\.slack: changes string to number \(declared by person/),
    ]);
  });

  it('rejects a redeclared reserved field', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/developer.yaml', '  owner?:        Person\n')))
      .toEqual(expect.arrayContaining([expect.stringMatching(/developer\.owner: reserved field name/)]));
  });

  it('rejects a prefix redeclared by a descendant', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/developer.yaml', 'prefix: P\n')))
      .toEqual(expect.arrayContaining([
        expect.stringMatching(/prefix P is declared by both (person and developer|developer and person)/),
        expect.stringMatching(/node developer: redeclares prefix inherited from person/),
      ]));
  });

  it('rejects an edge that widens its parent endpoints', () => {
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/COVERS.yaml'),
      fs.readFileSync(path.join(kg, 'edges/COVERS.yaml'), 'utf8').replace('to: Feature', 'to: Feature | Concept'))))
      .toEqual([expect.stringMatching(/COVERS\.to: Concept is not within EVIDENCES\.to: Feature/)]);
  });

  it('reports endpoint widening once, against the nearest violating ancestor', () => {
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/ASSERTS.yaml'),
      'type: ASSERTS\nextends: COVERS\nfrom: Scenario\nto: Concept\nderived_by: [annotation]\ndescription: Widens twice over.\n')))
      .toEqual(['edge ASSERTS.to: Concept is not within COVERS.to: Feature']);
  });

  it('quotes the ancestor declaration in a narrowing message without nesting quotes', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/developer.yaml', "  tier?:         \"'staff' | 'alien'\"\n")))
      .toEqual(["developer.tier: adds 'alien' not allowed by the ancestor (declared by person as 'staff' | 'contractor' | 'guest' = 'staff')"]);
  });

  it('rejects a duplicate YAML key', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/team.yaml', 'description: again\n')))
      .toEqual([expect.stringMatching(/YAML error: duplicated mapping key/)]);
  });

  it('rejects an unknown top-level key and a bad default type', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/team.yaml', 'colour: blue\n')))
      .toEqual([expect.stringMatching(/team: unknown key 'colour'/)]);
    expect(typeErrors((kg) => append(kg, 'nodes/team.yaml', "  size?:         \"'small' | 'large' = 'huge'\"\n")))
      .toEqual([expect.stringMatching(/team\.size: default 'huge' is not one of/)]);
  });

  it('rejects an edge key without annotation provenance and a same_type mismatch', () => {
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/COVERS.yaml'),
      fs.readFileSync(path.join(kg, 'edges/COVERS.yaml'), 'utf8').replace('derived_by: [annotation]', 'derived_by: [filesystem]'))))
      .toEqual([expect.stringMatching(/COVERS: has key 'covers' but 'annotation' is not in derived_by/)]);
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/PART_OF.yaml'),
      fs.readFileSync(path.join(kg, 'edges/PART_OF.yaml'), 'utf8').replace('to: Feature | Concept | Scenario', 'to: Feature'))))
      .toEqual([expect.stringMatching(/PART_OF: same_type but from/)]);
  });
});

describe('kg frontmatter', () => {
  it('splits a frontmatter block and reports where the body starts', () => {
    const fm = splitFrontmatter('---\nfeature: a\ntitle: T\n---\n\n# Heading\n');
    expect(fm.data).toEqual({feature: 'a', title: 'T'});
    expect(fm.bodyLine).toBe(5);
    expect(fm.body).toBe('\n# Heading\n');
  });

  it('treats a file without a leading --- as body only', () => {
    expect(splitFrontmatter('# Just a doc\n')).toMatchObject({data: null, bodyLine: 1});
  });

  it('rejects a duplicate key', () => {
    expect(splitFrontmatter('---\na: 1\na: 2\n---\n').error).toMatch(/duplicated mapping key/);
  });
});

describe('kg home documents (CONVENTIONS §5, §10)', () => {
  it('accepts the good monorepo: features, concepts, a developer and a team', () => {
    const repo = makeRepo();
    const {system, homes} = load(repo);
    expect(homes.errors).toEqual([]);
    const report = makeReport(system, homes);
    expect(report.homes).toEqual({feature: 5, concept: 2, developer: 1, team: 1});
    expect(report.unresolvedExternal).toBe(3);
    expect(report.stubs).toEqual(['platform', 'platform/caching']);
    expect(report.warnings.map((w) => w.message)).toEqual([
      expect.stringMatching(/no home for ~platform \(parent of ~platform\/caching\/invalidation\)/),
      expect.stringMatching(/no home for ~platform\/caching/),
    ]);
  });

  it('takes the name from title:, name: or the first heading', () => {
    const {homes} = load(makeRepo());
    const byId = new Map(homes.homes.map((h) => [h.id, h]));
    expect(byId.get('visualize/viewers/scatter-plot')!.name).toBe('Scatter plot');
    expect(byId.get('Team:core')!.name).toBe('Core team');
    expect(byId.get('C:column')!.name).toBe('Column');
  });

  it('resolves a reference through an alias', () => {
    const repo = makeRepo();
    const {homes} = load(repo);
    expect(homes.errors).toEqual([]);
    fs.writeFileSync(path.join(repo, 'public/help/visualize/viewers/scatter-plot.md'),
      fs.readFileSync(path.join(repo, 'public/help/visualize/viewers/scatter-plot.md'), 'utf8').replace('aliases: [visualize/scatterplot]', ''));
    expect(load(repo).homes.errors.map((e) => e.message))
      .toEqual([expect.stringMatching(/'visualize\/scatterplot' does not resolve to any home document/)]);
  });

  it('excludes Test Track scenarios, whose feature: key means something else', () => {
    const {homes} = load(makeRepo());
    expect(homes.homes.some((h) => h.file.includes('TestTrack'))).toBe(false);
    expect(homes.homes.some((h) => h.id === 'bio')).toBe(false);
  });

  it('rejects an unknown key', () => {
    expect(brokenErrors('unknown-key.md')).toEqual(["core/docs/broken.md:4: unknown key 'colour' for type feature"]);
  });

  it('requires owner on a second-level feature', () => {
    expect(brokenErrors('missing-owner.md')).toEqual(['core/docs/broken.md:2: ~govern/permissions is a second-level feature node and must set owner:']);
  });

  it('rejects a value outside its enum', () => {
    expect(brokenErrors('bad-enum.md')).toEqual(['core/docs/broken.md:4: status: "bogus" is not one of \'proposed\' | \'active\' | \'removed\'']);
  });

  it('rejects a reference that resolves nowhere, naming the id it looked for', () => {
    expect(brokenErrors('unresolved-ref.md')).toEqual(["core/docs/broken.md:4: concepts[0]: 'nosuch' does not resolve to any home document (as ~C:nosuch)"]);
  });

  it('rejects two homes claiming one id, with an error on each file naming the other', () => {
    expect(brokenErrors('duplicate-id.md')).toEqual([
      'core/docs/visualize.md:2: duplicate id ~visualize, also the home of core/docs/broken.md',
      'core/docs/broken.md:2: duplicate id ~visualize, also the home of core/docs/visualize.md',
    ]);
  });

  it('finds a home anywhere in the repos, not only under the docs folders', () => {
    const repo = makeRepo();
    const doc = '---\nfeature: visualize/viewers/box-plot\n---\n# Box plot\n';
    for (const file of ['public/packages/Chem/docs/box-plot.md', 'core/tools/grok-core/tool-doc.md', 'infra/docs/hosts.md', 'landing/pages/pricing.mdx'])
      write(repo, file, doc.replace('box-plot', path.basename(file).replace(/\.mdx?$/, '')));
    for (const file of ['public/packages/Chem/node_modules/x/README.md', 'core/.git/README.md', 'public/help/dist/x.md',
      'core/tools/grok-core/fixtures/reports/1.md', 'public/tools/bin/__tests__/fixtures/kg/good/core/docs/x.md',
      'public/playwright-public/Bio/analyze.md', 'public/packages/UsageAnalysis/files/_template-ui.md'])
      write(repo, file, doc.replace('box-plot', 'hidden'));
    const {homes} = load(repo);
    expect(homes.errors).toEqual([]);
    expect(homes.homes.map((h) => h.id).filter((id) => id.startsWith('visualize/viewers/')).sort())
      .toEqual(['box-plot', 'hosts', 'old-scatter', 'pricing', 'scatter-plot', 'tool-doc'].map((s) => `visualize/viewers/${s}`));
    expect(homes.homes.some((h) => h.id.endsWith('hidden'))).toBe(false);
  });

  it('rejects a backslash in any Path: code items, path maps and body citations', () => {
    const repo = makeRepo();
    write(repo, 'core/docs/broken.md', '---\nfeature: govern/permissions\nowner: askalkin\ncode:\n  - core\\client\\d4\n  - {path: "public\\\\js-api", role: api}\n---\n# Permissions\n\nSee `core\\client\\d4\\lib\\scatter.dart`.\n');
    expect(load(repo).homes.errors.map((e) => e.message)).toEqual([
      "code[0]: path 'core\\client\\d4' contains a backslash; use forward slashes",
      "code[1]: path 'public\\js-api' contains a backslash; use forward slashes",
      "cited path 'core\\client\\d4\\lib\\scatter.dart' contains a backslash; use forward slashes",
    ]);
  });

  it('accepts title: as the name source outside public/help too', () => {
    const repo = makeRepo();
    write(repo, 'core/docs/broken.md', '---\nfeature: govern/permissions\nowner: askalkin\ntitle: Permissions\n---\nNo heading here.\n');
    const {homes} = load(repo);
    expect(homes.errors).toEqual([]);
    expect(homes.homes.find((h) => h.id === 'govern/permissions')!.name).toBe('Permissions');
  });

  it('reports a name: of the wrong kind once, without the "no name" hint', () => {
    const repo = makeRepo();
    write(repo, 'core/docs/broken.md', '---\nfeature: govern/permissions\nowner: askalkin\nname: 5\n---\nNo heading here.\n');
    expect(load(repo).homes.errors.map((e) => e.message)).toEqual(['name: expected a string, got 5']);
  });

  it('tells the author to quote an unquoted Version without inventing a number', () => {
    const repo = makeRepo();
    write(repo, 'core/docs/broken.md', '---\nfeature: govern/permissions\nowner: askalkin\nsuperseded_by: [{to: visualize, release: 1.20}]\n---\n# Permissions\n');
    expect(load(repo).homes.errors.map((e) => e.message))
      .toEqual(['superseded_by[0].release: a Version must be a quoted string; YAML read the unquoted value as a number, quote the value']);
  });

  it('rejects a cited path that does not exist, with its body line', () => {
    expect(brokenErrors('cited-path-missing.md')).toEqual(["core/docs/broken.md:7: cited path 'core/server/nope.dart' does not exist"]);
  });

  it('validates edge properties on a map item against the edge type', () => {
    expect(brokenErrors('bad-edge-property.md')).toEqual([
      'core/docs/broken.md:4: concepts[0].role: "bogus" is not one of \'central\' | \'supporting\'',
      "core/docs/broken.md:4: concepts[1]: USES_CONCEPT has no property 'weight' (role)",
    ]);
  });

  it('rejects an authored part_of', () => {
    expect(brokenErrors('part-of.md')).toEqual(['core/docs/broken.md:4: PART_OF is derived from the id path and never authored; remove part_of:']);
  });

  it('rejects an edge key on a home whose type is not the edge endpoint', () => {
    expect(brokenErrors('wrong-endpoint.md')).toEqual(['core/docs/broken.md:3: tickets: TRACKED_IN.from must be Feature; this home is a concept']);
  });

  it('rejects type: that is not a subtype of the prefix type, an abstract type, and an extracted type', () => {
    expect(brokenErrors('bad-subtype.md')).toEqual(["core/docs/broken.md:3: type 'customer' is not a subtype of person, the type of prefix P"]);
    expect(brokenErrors('abstract-type.md')).toEqual(["core/docs/broken.md:3: type 'actor' is not a subtype of person, the type of prefix P"]);
    expect(brokenErrors('extracted-type.md')).toEqual(["core/docs/broken.md:2: type 'release' is extracted, not authored; it cannot have a home document"]);
  });

  it('checks the edges: escape hatch for authorable types, endpoints and existence', () => {
    expect(brokenErrors('edges-escape-hatch.md')).toEqual([
      'core/docs/broken.md:4: edges[0]: PART_OF is never authored (derived_by lacks annotation)',
      'core/docs/broken.md:4: edges[1]: COVERS.from must be Scenario; this home is a feature',
      'core/docs/broken.md:4: edges[2]: unknown edge type "NOPE"',
    ]);
  });

  it('rejects a duplicate frontmatter key', () => {
    expect(brokenErrors('duplicate-yaml-key.md')).toEqual(['core/docs/broken.md:4: YAML error: duplicated mapping key']);
  });
});

describe('kg gen (CONVENTIONS §11.2)', () => {
  const system = loadTypeSystem(path.join(fixtures, 'good', KG_DIR));

  it('writes kg.d.ts as an interface hierarchy over Ref<T>', () => {
    const dts = generateDts(system);
    expect(dts).toContain('type Ref<T> = string;');
    expect(dts).toContain('export interface Developer extends Person {\n  company: Ref<Team>;\n  bitbucket: string;\n  areas?: Ref<Feature>[];\n}');
    expect(dts).toContain("  owner?: Ref<Person | Team>;");
    expect(dts).toContain('export interface CoversEdge {\n  from: Ref<Scenario>;\n  to: Ref<Feature>;\n  derived_by: Provenance;\n  confidence: number;\n  evidence?: Path[];\n  strength?: \'weak\' | \'normal\' | \'strong\';\n  level?: \'exercised\' | \'asserted\';\n}');
    expect(dts).not.toContain('interface EvidencesEdge');
    expect(dts).toContain("export type NodeTypeName = 'actor' | 'artifact'");
    expect(generateDts(system)).toBe(dts);
  });

  it('splices the glossary tables and keeps the Concepts section byte for byte', () => {
    const before = fs.readFileSync(path.join(fixtures, 'good', KG_DIR, 'GLOSSARY.md'), 'utf8');
    const after = spliceGlossary(before, system);
    expect(after.slice(after.indexOf('## Concepts'))).toBe(before.slice(before.indexOf('## Concepts')));
    expect(after.slice(0, after.indexOf('## Prefixes'))).toBe(before.slice(0, before.indexOf('## Prefixes')));
    expect(after).not.toContain('stale');
    expect(after).toContain('| `~C:a/b` | concept | yes | yes | core/docs/concepts/<name>.md |');
    expect(after).toContain('| `~Rel:name` | release | no | no | extracted |');
    expect(after).toContain('| developer | actor | person |  | A person who commits to the platform. |');
    expect(after).toContain('| `PART_OF` | Feature \\| Concept \\| Scenario → Feature \\| Concept \\| Scenario |  |  | filesystem |');
    expect(after).toContain('| `SUPERSEDES` | Feature → Feature |  | `superseded_by:` (on target) | annotation |');
  });

  it('renders FEATURES.md as an indented tree with stubs for missing parents', () => {
    const repo = makeRepo();
    const {homes} = load(repo);
    const md = generateFeatures(system, homes);
    expect(md.split('\n').slice(4)).toEqual([
      '## Features',
      '',
      '- `~platform` *(no home yet)*',
      '  - `~platform/caching` *(no home yet)*',
      '    - `~platform/caching/invalidation` — Cache invalidation ([home](CACHING.md))',
      '- `~visualize` — Visualize ([home](visualize.md))',
      '  - `~visualize/viewers` — Viewers ([home](viewers/README.md))',
      '    - `~visualize/viewers/old-scatter` — Old scatter ([home](viewers/old-scatter.md))',
      '    - `~visualize/viewers/scatter-plot` — Scatter plot ([home](../../public/help/visualize/viewers/scatter-plot.md))',
      '',
      '## Concepts',
      '',
      '- `~C:column` — Column ([home](concepts/column.md))',
      '- `~C:dataframe` — Dataframe ([home](concepts/dataframe.md))',
      '',
    ]);
  });

  it('writes a note instead of a tree when there are no homes', () => {
    const md = generateFeatures(system, {homes: [], stubs: [], errors: [], warnings: [], unresolvedExternal: 0, scanned: 0});
    expect(md).toContain('No home documents yet');
  });

  it('reports a GLOSSARY.md without the expected headings as an error instead of throwing', () => {
    const repo = makeRepo();
    const kgRoot = path.join(repo, KG_DIR);
    const glossary = path.join(kgRoot, 'GLOSSARY.md');
    fs.writeFileSync(glossary, fs.readFileSync(glossary, 'utf8').replace('## Concepts', '## Terms'));
    const {outputs, errors} = generate(system, kgRoot, repo, load(repo).homes);
    expect(errors).toEqual([{file: glossary, message: expect.stringMatching(/^heading '## Concepts' not found; the generated tables are spliced between/)}]);
    expect(outputs.map((o) => path.basename(o.file))).toEqual(['kg.d.ts', 'FEATURES.md']);
  });

  it('gen --check detects drift and writes nothing', () => {
    const repo = makeRepo();
    const kgRoot = path.join(repo, KG_DIR);
    const {homes} = load(repo);
    const {outputs} = generate(system, kgRoot, repo, homes);
    expect(writeOutputs(outputs, true).stale.map((f) => path.basename(f))).toEqual(['kg.d.ts', 'GLOSSARY.md', 'FEATURES.md']);
    expect(fs.existsSync(path.join(kgRoot, 'kg.d.ts'))).toBe(false);
    expect(writeOutputs(outputs, false).written).toHaveLength(3);
    expect(writeOutputs(outputs, true).stale).toEqual([]);
    fs.appendFileSync(path.join(kgRoot, 'kg.d.ts'), '// drift\n');
    expect(writeOutputs(outputs, true).stale.map((f) => path.basename(f))).toEqual(['kg.d.ts']);
  });
});

describe('grok kg command', () => {
  it('check --output json prints the report and nothing else', async () => {
    const repo = makeRepo();
    const {ok, out, exitCode} = await run({_: ['kg', 'check'], kg: path.join(repo, KG_DIR), output: 'json'});
    expect(ok).toBe(true);
    expect(exitCode).toBeUndefined();
    expect(out).toHaveLength(1);
    const report = JSON.parse(out[0]);
    expect(report.errors).toEqual([]);
    expect(report.types).toEqual({nodes: 17, edges: 8, prefixes: 6});
    expect(report.homes.feature).toBe(5);
    expect(report.scanned).toBeGreaterThan(9);
  });

  it('--types-only skips the home documents', async () => {
    const repo = makeRepo();
    const {out} = await run({_: ['kg', 'check'], kg: path.join(repo, KG_DIR), output: 'json', 'types-only': true});
    expect(JSON.parse(out[0]).homes).toEqual({});
  });

  it('names an unknown verb and returns false so grok.js prints the usage', async () => {
    const {ok, err} = await run({_: ['kg', 'frobnicate']});
    expect(ok).toBe(false);
    expect(err).toEqual(["unknown verb 'frobnicate'"]);
  });

  it('refuses planned verbs, a bad --output, a stray argument and gen --types-only with exit 1', async () => {
    const repo = makeRepo();
    const kgRoot = path.join(repo, KG_DIR);
    for (const [argv, message] of [
      [{_: ['kg', 'build']}, /grok kg build is not implemented yet/],
      [{_: ['kg', 'check'], kg: kgRoot, output: 'csv'}, /--output must be table or json, got 'csv'/],
      [{_: ['kg', 'check', 'extra'], kg: kgRoot}, /unexpected argument 'extra'/],
      [{_: ['kg', 'gen'], kg: kgRoot, check: true, 'types-only': true}, /--types-only cannot be combined with gen/],
      [{_: ['kg', 'check'], kg: 'C:\\nowhere\\kg'}, /^C:\/nowhere\/kg: no schema\.yaml$/],
    ] as const) {
      const {ok, out, err, exitCode} = await run(argv);
      expect(ok).toBe(true);
      expect(exitCode).toBe(1);
      expect(out).toEqual([]);
      expect(err).toEqual([expect.stringMatching(message)]);
    }
  });
});
