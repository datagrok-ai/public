/// `grok kg check` / `grok kg gen` against the mini monorepo under fixtures/kg/good and the
/// deliberately broken home documents under fixtures/kg/broken (conventions.md §5, §7, §10).
import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {Project} from 'ts-morph';
import {parseMember, loadTypeSystem, graphLabel, TypeSystem} from '../utils/kg/types';
import {splitFrontmatter} from '../utils/kg/frontmatter';
import {extractCitations, proseLines} from '../utils/kg/citations';
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

describe('kg member syntax (conventions.md §7.1)', () => {
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
    expect(parseMember('t?', "'a' | 'b' = 'c'").error).toMatch(/default "c": "c" is not one of 'a' \| 'b'/);
    expect(parseMember('n?', "number = 'x'").error).toMatch(/expected a finite number/);
  });

  it('rejects a bad member name and a non-string spec', () => {
    expect(parseMember('Bad-Name', 'string').error).toMatch(/bad member name/);
    expect(parseMember('n', 5).error).toMatch(/must be a string/);
  });
});

describe('kg type files (conventions.md §7.3, §7.4, schema.yaml constraints)', () => {
  it('validates the fixture type system clean and merges members down the chain', () => {
    const system = loadTypeSystem(path.join(fixtures, 'good', KG_DIR));
    expect(system.errors).toEqual([]);
    const developer = system.nodes.get('developer')!;
    expect(developer.chain).toEqual(['developer', 'person', 'actor', 'node']);
    expect(developer.prefix).toBe('P');
    expect(developer.authored).toBe(true);
    expect(developer.members.company).toMatchObject({nullable: false, refs: ['team']});
    expect(developer.members.email).toMatchObject({nullable: false});
    expect(system.edges.get('covers')!.properties).toHaveProperty('strength');
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
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/covers.yaml'),
      fs.readFileSync(path.join(kg, 'edges/covers.yaml'), 'utf8').replace('to: Feature', 'to: Feature | Concept'))))
      .toEqual([expect.stringMatching(/covers\.to: Concept is not within evidences\.to: Feature/)]);
  });

  it('reports endpoint widening once, against the nearest violating ancestor', () => {
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/asserts.yaml'),
      'type: asserts\nextends: covers\nfrom: Scenario\nto: Concept\nderived_by: [annotation]\ndescription: Widens twice over.\n')))
      .toEqual(['edge asserts.to: Concept is not within covers.to: Feature']);
  });

  it('quotes the ancestor declaration in a narrowing message without nesting quotes', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/developer.yaml', "  tier?:         \"'staff' | 'alien'\"\n")))
      .toEqual(["developer.tier: adds 'alien' not allowed by the ancestor (declared by person as 'staff' | 'contractor' | 'guest' = 'staff')"]);
  });

  it('rejects an upper-snake type name and points at the lower-dash-case spelling', () => {
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/DEPENDS_ON.yaml'),
      'type: DEPENDS_ON\nfrom: Feature\nto: Feature\nderived_by: [annotation]\ndescription: Old spelling.\n')))
      .toEqual([expect.stringMatching(/^edge type name 'DEPENDS_ON' is upper-snake; type names are lower-dash-case, write 'depends-on' \(the graph label DEPENDS_ON is derived from it\)$/)]);
  });

  it('renders the graph label of a type name as upper-snake', () => {
    expect(graphLabel('part-of')).toBe('PART_OF');
    expect(graphLabel('uses-concept')).toBe('USES_CONCEPT');
    expect(graphLabel('covers')).toBe('COVERS');
  });

  it('rejects a duplicate YAML key', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/team.yaml', 'description: again\n')))
      .toEqual([expect.stringMatching(/YAML error: duplicated mapping key/)]);
  });

  it('rejects an unknown top-level key and a bad default type', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/team.yaml', 'colour: blue\n')))
      .toEqual([expect.stringMatching(/team: unknown key 'colour'/)]);
    expect(typeErrors((kg) => append(kg, 'nodes/team.yaml', "  size?:         \"'small' | 'large' = 'huge'\"\n")))
      .toEqual([expect.stringMatching(/team\.size: default "huge": "huge" is not one of/)]);
  });

  it('rejects an edge key without annotation provenance and a same_type mismatch', () => {
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/covers.yaml'),
      fs.readFileSync(path.join(kg, 'edges/covers.yaml'), 'utf8').replace('derived_by: [annotation]', 'derived_by: [filesystem]'))))
      .toEqual([expect.stringMatching(/covers: has key 'covers' but 'annotation' is not in derived_by/)]);
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/part-of.yaml'),
      fs.readFileSync(path.join(kg, 'edges/part-of.yaml'), 'utf8').replace('to: Feature | Concept | Scenario', 'to: Feature'))))
      .toEqual([expect.stringMatching(/part-of: same_type but from/)]);
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

describe('kg home documents (conventions.md §5, §10)', () => {
  it('accepts the good monorepo: features, concepts, a developer and a team', () => {
    const repo = makeRepo();
    const {system, homes} = load(repo);
    expect(homes.errors).toEqual([]);
    const report = makeReport(system, homes);
    expect(report.homes).toEqual({feature: 6, concept: 2, developer: 1, team: 1, scenario: 1});
    expect(report.annotatedPages).toBe(1);
    expect(report.unresolvedExternal).toHaveLength(4);
    expect(report.unresolvedExternal).toContainEqual({source: 'core/docs/viewers/README.md', key: 'tickets[0]', value: 'GROK-20863', expectedTypes: ['Ticket']});
    expect(report.stubs).toEqual(['platform']);
    expect(report.warnings.map((w) => w.message)).toEqual([expect.stringMatching(/no home for ~platform \(parent of ~platform\/caching\)/)]);
    expect(report.codes).toEqual({stub: 1});
    expect(report.citations).toEqual({doc: 1, code: 3});
  });

  it('takes the name from title:, name: or the first heading', () => {
    const {homes} = load(makeRepo());
    const byId = new Map(homes.homes.map((h) => [h.id, h]));
    expect(byId.get('visualize/viewers/scatter-plot')!.name).toBe('Scatter plot');
    expect(byId.get('Team:core')!.name).toBe('Core team');
    expect(byId.get('C:column')!.name).toBe('Column');
  });

  it('reads a YAML record as a home: same keys, edge-key map items, prose in description:', () => {
    const {homes} = load(makeRepo());
    const dataframe = homes.homes.find((h) => h.id === 'C:dataframe')!;
    expect(homes.errors).toEqual([]);
    expect(dataframe).toMatchObject({yaml: true, file: 'core/docs/knowledge-graph/concepts/dataframe.yaml', name: 'DataFrame', line: 1});
    expect(homes.homes.filter((h) => h.yaml).map((h) => h.id).sort()).toEqual(['C:column', 'C:dataframe', 'P:askalkin', 'Team:core']);
  });

  it('requires name: on a YAML home and reports the key line of a bad map-item property', () => {
    const repo = makeRepo();
    write(repo, 'core/docs/knowledge-graph/concepts/entity.yaml', 'id: C:entity\ndescription: Any server-side object.\n');
    write(repo, 'core/docs/knowledge-graph/concepts/query.yaml', 'id: C:query\nname: Query\ndefined_by:\n  - {to: decl:core/client/d4/lib/scatter.dart#Query, language: cobol}\n');
    expect(load(repo).homes.errors.map((e) => `${e.file}:${e.line}: ${e.message}`)).toEqual([
      'core/docs/knowledge-graph/concepts/entity.yaml:1: no name: a YAML home must set name:',
      'core/docs/knowledge-graph/concepts/query.yaml:3: defined_by[0].language: "cobol" is not one of \'dart\' | \'ts\'',
    ]);
  });

  it('warns about a YAML file in the knowledge-graph folder that is not a home, and skips nodes/, edges/ and schema.yaml', () => {
    const repo = makeRepo();
    write(repo, 'core/docs/knowledge-graph/notes.yaml', 'todo: [migrate scenarios]\n');
    const {homes} = load(repo);
    expect(homes.errors).toEqual([]);
    expect(homes.warnings.map((w) => `${w.file}: ${w.message}`)).toContain(
      'core/docs/knowledge-graph/notes.yaml: stray YAML file in the knowledge-graph folder: no id: or feature: key, so not a home');
    expect(homes.warnings.some((w) => /nodes\/|edges\/|schema\.yaml/.test(w.file))).toBe(false);
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

  it('skips legacy Test Track files, whose feature: key means the area, but reads migrated id: TS: ones', () => {
    const {homes} = load(makeRepo());
    expect(homes.homes.filter((h) => h.file.includes('TestTrack')).map((h) => h.id)).toEqual(['TS:scatter-plot-ui']);
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
      "core/docs/broken.md:4: concepts[1]: uses-concept has no property 'weight' (role)",
    ]);
  });

  it('rejects an authored part_of', () => {
    expect(brokenErrors('part-of.md')).toEqual(['core/docs/broken.md:4: part-of is derived from the id path and never authored; remove part_of:']);
  });

  it('rejects an edge key on a home whose type is not the edge endpoint', () => {
    expect(brokenErrors('wrong-endpoint.md')).toEqual(['core/docs/broken.md:3: tickets: tracked-in.from must be Feature; this home is a concept']);
  });

  it('rejects type: that is not a subtype of the prefix type, an abstract type, and an extracted type', () => {
    expect(brokenErrors('bad-subtype.md')).toEqual(["core/docs/broken.md:3: type 'customer' is not a subtype of person, the type of prefix P"]);
    expect(brokenErrors('abstract-type.md')).toEqual(["core/docs/broken.md:3: type 'actor' is not a subtype of person, the type of prefix P"]);
    expect(brokenErrors('extracted-type.md')).toEqual(["core/docs/broken.md:2: type 'release' is extracted, not authored; it cannot have a home document"]);
  });

  it('checks the edges: escape hatch for authorable types, endpoints and existence', () => {
    expect(brokenErrors('edges-escape-hatch.md')).toEqual([
      'core/docs/broken.md:4: edges[0]: part-of is never authored (derived_by lacks annotation)',
      'core/docs/broken.md:4: edges[1]: covers.from must be Scenario; this home is a feature',
      'core/docs/broken.md:4: edges[2]: unknown edge type "nope"',
      "core/docs/broken.md:4: edges[3]: edge types are lower-dash-case: write 'uses-concept', not 'USES_CONCEPT'",
    ]);
  });

  it('rejects a duplicate frontmatter key', () => {
    expect(brokenErrors('duplicate-yaml-key.md')).toEqual(['core/docs/broken.md:4: YAML error: duplicated mapping key']);
  });
});

describe('kg gen (conventions.md §11.2)', () => {
  const system = loadTypeSystem(path.join(fixtures, 'good', KG_DIR));

  it('writes kg.d.ts as an interface hierarchy over Ref<T>', () => {
    const dts = generateDts(system);
    expect(dts).toContain('type Ref<T> = string;');
    expect(dts).toContain('export interface Developer extends Person {\n  company: Ref<Team>;\n  bitbucket: string;\n  areas?: Ref<Feature>[];\n}');
    expect(dts).toContain("  owner?: Ref<Person | Team>;");
    expect(dts).toContain('export interface CoversEdge {\n  from: Ref<Scenario>;\n  to: Ref<Feature>;\n  derived_by: Provenance;\n  confidence: number;\n  evidence?: Path[];\n  batch: string;\n  strength?: \'weak\' | \'normal\' | \'strong\';\n  level?: \'exercised\' | \'asserted\';\n}');
    expect(dts).not.toContain('interface EvidencesEdge');
    expect(dts).toContain("export type NodeTypeName = 'actor' | 'artifact'");
    expect(dts).toContain('/** The hierarchy edge, derived from the id path and never authored. */\nexport interface PartOfEdge {');
    expect(dts).toContain('export interface UsesConceptEdge {');
    expect(dts).toContain("export type EdgeTypeName = 'covers' | 'defines-concept' | 'documents' | 'is-implemented-in' | 'mentions' | 'part-of' | 'supersedes' | 'tracked-in' | 'uses-concept';");
    expect(dts).toContain("export type RefPredicate = 'areas' | 'company' | 'lead' | 'owner';");
    expect(dts).toContain('  aliases?: string[];\n  source_layer: \'public\' | \'core\' | \'infra\' | \'process\' | \'synthetic\';\n  home?: Path;\n  provenance: Provenance;\n  batch: string;\n}');
    expect(dts).toContain('  to: Ref<Feature>;\n  derived_by: Provenance;\n  confidence: number;\n  evidence?: Path[];\n  batch: string;\n');
    expect(generateDts(system)).toBe(dts);
  });

  it('splices the glossary tables; without homes the Concepts section is kept byte for byte', () => {
    const before = fs.readFileSync(path.join(fixtures, 'good', KG_DIR, 'glossary.md'), 'utf8');
    const after = spliceGlossary(before, system);
    expect(after.slice(after.indexOf('## Concepts'))).toBe(before.slice(before.indexOf('## Concepts')));
    expect(after.slice(0, after.indexOf('## Prefixes'))).toBe(before.slice(0, before.indexOf('## Prefixes')));
    expect(after).not.toContain('stale');
    expect(after).toContain('| `~C:a/b` | concept | yes | yes | core/docs/knowledge-graph/concepts/<name>.yaml |');
    expect(after).toContain('| `~Rel:name` | release | no | no | extracted |');
    expect(after).toContain('| developer | actor | person |  | A person who commits to the platform. |');
    expect(after).toContain('| Edge | Label | From → To | Extends | Key | Derived by | One line |');
    expect(after).toContain('| `part-of` | PART_OF | Feature \\| Concept \\| Scenario → Feature \\| Concept \\| Scenario |  |  | filesystem |');
    expect(after).toContain('| `supersedes` | SUPERSEDES | Feature → Feature |  | `superseded_by:` (on target) | annotation |');
    expect(after).toContain('| `covers` | COVERS | Scenario → Feature | evidences | `covers:` | annotation |');
  });

  it('generates the Concepts table from the concept homes, keeping the intro lines above it', () => {
    const before = fs.readFileSync(path.join(fixtures, 'good', KG_DIR, 'glossary.md'), 'utf8');
    const after = spliceGlossary(before, system, load(makeRepo()).homes);
    expect(after.slice(after.indexOf('## Concepts'))).toBe([
      '## Concepts',
      '',
      'Hand-written concepts table, byte for byte.',
      '',
      '| Concept | Name | Area | One line | Defined by |',
      '|---|---|---|---|---|',
      '| `column` | Column |  | A typed vector inside a dataframe. |  |',
      '| `dataframe` | DataFrame | data | An in-memory columnar table. | `DataFrame` |',
      '',
    ].join('\n'));
    expect(after).not.toContain('| A table |');
  });

  it('renders feature-tree.md as an indented tree with links relative to the knowledge-graph folder', () => {
    const repo = makeRepo();
    const {homes} = load(repo);
    const md = generateFeatures(system, homes);
    expect(md.split('\n').slice(4)).toEqual([
      '## Features',
      '',
      '- `~platform` *(no home yet)*',
      '  - `~platform/caching` — Caching ([home](../CACHING.md))',
      '    - `~platform/caching/invalidation` — Cache invalidation ([home](../caching/invalidation.md))',
      '- `~visualize` — Visualize ([home](../visualize.md))',
      '  - `~visualize/viewers` — Viewers ([home](../viewers/README.md))',
      '    - `~visualize/viewers/old-scatter` — Old scatter ([home](../viewers/old-scatter.md))',
      '    - `~visualize/viewers/scatter-plot` — Scatter plot ([home](../../../public/help/visualize/viewers/scatter-plot.md))',
      '',
      '## Concepts',
      '',
      '- `~C:column` — Column ([home](concepts/column.yaml))',
      '- `~C:dataframe` — DataFrame ([home](concepts/dataframe.yaml))',
      '',
      '## Scenarios',
      '',
      '- `~TS:scatter-plot-ui` — Scatter plot manual checks ([home](../../../public/packages/UsageAnalysis/files/TestTrack/Viewers/scatter-plot-ui.md))',
      '',
    ]);
  });

  it('writes a note instead of a tree when there are no homes', () => {
    const md = generateFeatures(system, {homes: [], stubs: [], errors: [], warnings: [], unresolvedExternal: [], scanned: 0, annotatedPages: 0, citations: {doc: 0, code: 0}});
    expect(md).toContain('No home documents yet');
  });

  it('reports a glossary.md without the expected headings as an error instead of throwing', () => {
    const repo = makeRepo();
    const kgRoot = path.join(repo, KG_DIR);
    const glossary = path.join(kgRoot, 'glossary.md');
    fs.writeFileSync(glossary, fs.readFileSync(glossary, 'utf8').replace('## Concepts', '## Terms'));
    const {outputs, errors} = generate(system, kgRoot, repo, load(repo).homes);
    expect(errors).toEqual([{file: glossary, code: 'glossary', message: expect.stringMatching(/^heading '## Concepts' not found; the generated tables are spliced between/)}]);
    expect(outputs.map((o) => path.basename(o.file))).toEqual(['kg.d.ts', 'feature-tree.md']);
  });

  it('gen --check detects drift and writes nothing', () => {
    const repo = makeRepo();
    const kgRoot = path.join(repo, KG_DIR);
    const {homes} = load(repo);
    const {outputs} = generate(system, kgRoot, repo, homes);
    expect(writeOutputs(outputs, true).stale.map((f) => path.basename(f))).toEqual(['kg.d.ts', 'glossary.md', 'feature-tree.md']);
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
    expect(report.types).toEqual({nodes: 19, edges: 11, prefixes: 6});
    expect(report.homes).toMatchObject({feature: 6, scenario: 1});
    expect(report.annotatedPages).toBe(1);
    expect(Array.isArray(report.unresolvedExternal)).toBe(true);
    expect(report.codes).toEqual({stub: 1});
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

  it('refuses a nameless report, a bad --output, a stray argument and gen --types-only with exit 1', async () => {
    const repo = makeRepo();
    const kgRoot = path.join(repo, KG_DIR);
    for (const [argv, message] of [
      [{_: ['kg', 'report'], kg: kgRoot}, /^grok kg report needs a report name: orphans, stale, coverage, proposed, diff$/],
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

/** Errors of a copy of the good repo after a home document was dropped in at [file]. */
function homeErrors(file: string, text: string, mutate?: (repo: string) => void): string[] {
  const repo = makeRepo();
  write(repo, file, text);
  mutate?.(repo);
  return load(repo).homes.errors.map((e) => `${e.code} ${e.file}:${e.line}: ${e.message}`);
}

const PERMISSIONS = '---\nfeature: govern/permissions\nowner: askalkin\n';

describe('reference resolution: parse, expand, type-check, look up (review 2 #1)', () => {
  it('rejects ids whose kind cannot be of the expected type instead of calling them external', () => {
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}owner: GROK-123\n---\n# P\n`.replace('owner: askalkin\n', '')))
      .toEqual(["unresolved-ref core/docs/broken.md:3: owner: 'GROK-123' is a Ticket; expected Person | Team"]);
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}concepts: [Rel:1.28]\n---\n# P\n`))
      .toEqual(["unresolved-ref core/docs/broken.md:4: concepts[0]: 'Rel:1.28' is a release (prefix Rel); expected Concept"]);
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}concepts: [Xyz:foo]\n---\n# P\n`))
      .toEqual(["unresolved-ref core/docs/broken.md:4: concepts[0]: 'Xyz:foo': unknown id prefix 'Xyz'"]);
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}tickets: [some-thing]\n---\n# P\n`))
      .toEqual([expect.stringMatching(/tickets\[0\]: 'some-thing' cannot be a Ticket: an extracted id carries its scheme or tracker key/)]);
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}concepts: [Not Valid]\n---\n# P\n`))
      .toEqual([expect.stringMatching(/'Not Valid' is not a valid id/)]);
  });

  it('expands an abstract reference type to its concrete authored descendants and their prefixes', () => {
    const repo = makeRepo();
    fs.appendFileSync(path.join(repo, KG_DIR, 'nodes/team.yaml'), '  sponsor?:      Actor\n');
    write(repo, 'core/docs/knowledge-graph/internal/teams/core.yaml', 'id: Team:core\nname: Core team\nsponsor: askalkin\n');
    expect(load(repo).homes.errors).toEqual([]);
    write(repo, 'core/docs/knowledge-graph/internal/teams/core.yaml', 'id: Team:core\nname: Core team\nsponsor: nobody\n');
    expect(load(repo).homes.errors.map((e) => e.message))
      .toEqual(["sponsor: 'nobody' does not resolve to any home document (as ~Cust:nobody or ~P:nobody or ~Team:nobody)"]);
  });

  it('checks the file part of a decl: reference now and defers only the symbol', () => {
    expect(homeErrors('core/docs/knowledge-graph/concepts/entity.yaml', 'id: C:entity\nname: Entity\ndefined_by: [decl:core/nope.dart#Entity]\n'))
      .toEqual(["unresolved-ref core/docs/knowledge-graph/concepts/entity.yaml:3: defined_by[0]: 'decl:core/nope.dart#Entity': the declaration's path 'core/nope.dart' does not exist"]);
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}concepts: [what:ever]\n---\n# P\n`))
      .toEqual(["unresolved-ref core/docs/broken.md:4: concepts[0]: 'what:ever': unknown id scheme 'what:'"]);
  });
});

describe('annotated pages and migrated scenarios (review 2 #2)', () => {
  it('validates a page that only carries documents: and counts it', () => {
    const repo = makeRepo();
    write(repo, 'public/help/visualize/viewers/more-tips.md', '---\ntitle: More\ndocuments:\n  - nosuch\n  - {to: visualize/viewers/scatter-plot, audience: robot}\ncovers: [visualize]\n---\nText.\n');
    const {homes} = load(repo);
    expect(homes.annotatedPages).toBe(2);
    expect(homes.errors.map((e) => `${e.file}:${e.line}: ${e.message}`)).toEqual([
      "public/help/visualize/viewers/more-tips.md:3: documents[0]: 'nosuch' does not resolve to any home document (as ~nosuch)",
      'public/help/visualize/viewers/more-tips.md:3: documents[1].audience: "robot" is not one of \'user\' | \'developer\'',
      'public/help/visualize/viewers/more-tips.md:6: covers: covers.from must be Scenario; this page is a doc-page',
    ]);
  });

  it('reads a migrated Test Track scenario as a home, derives path and manual_only, rejects an authored path', () => {
    const {homes} = load(makeRepo());
    const scenario = homes.homes.find((h) => h.id === 'TS:scatter-plot-ui')!;
    expect(scenario.data.path).toBe('public/packages/UsageAnalysis/files/TestTrack/Viewers/scatter-plot-ui.md');
    expect(scenario.data.manual_only).toBe(true);
    expect(homes.homes.some((h) => h.id === 'bio')).toBe(false);
    expect(homeErrors('public/packages/UsageAnalysis/files/TestTrack/Viewers/other.md', '---\nid: TS:other\npath: elsewhere.md\n---\n# Other\n'))
      .toEqual(['authored-path public/packages/UsageAnalysis/files/TestTrack/Viewers/other.md:3: path: derived from the home file (public/packages/UsageAnalysis/files/TestTrack/Viewers/other.md) and never authored; remove it']);
  });
});

describe('citations as records (review 2 #4)', () => {
  it('extracts backticks, inline links relative to the document, reference-style links and in-page anchors', () => {
    const body = [
      'See `core/client/d4/lib/scatter.dart:12` and [tips](scatter-plot-tips.md#tips) and [root](/core/docs/CACHING.md).',
      'Also [spaced](<../viewers/my file.md>) and [ref][spec] and [own](#usage).',
      '',
      '[spec]: ../../../../core/docs/SPEC.md',
      '',
      '````md',
      '```',
      'core/inside/fence.md',
      '```',
      '````',
      '[out](../../../../../etc/passwd)',
    ].join('\n');
    const cites = extractCitations('public/help/visualize/viewers/scatter-plot.md', body, 10);
    expect(cites.map((c) => [c.kind, c.resolved, c.anchor ?? '', c.target, c.line])).toEqual([
      ['backtick', 'core/client/d4/lib/scatter.dart', '', 'code', 10],
      ['link', 'public/help/visualize/viewers/scatter-plot-tips.md', 'tips', 'doc', 10],
      ['link', 'core/docs/CACHING.md', '', 'doc', 10],
      ['link', 'public/help/visualize/viewers/my file.md', '', 'doc', 11],
      ['link', 'public/help/visualize/viewers/scatter-plot.md', 'usage', 'doc', 11],
      ['reference-link', 'core/docs/SPEC.md', '', 'doc', 13],
      ['link', null, '', 'code', 20],
    ]);
    expect(proseLines('````\n```\ninside\n```\n````\nout').map((l) => l.text)).toEqual(['out']);
  });

  it('reports doc links, code citations, anchors and repository escapes separately', () => {
    const errors = homeErrors('core/docs/broken.md', `${PERMISSIONS}---\n# P\n\n[a](../nope.md) [b](viewers/README.md#no-such) \`core/nope.dart\` [c](../../../x.md) [d](#missing)\n`);
    expect(errors).toEqual([
      "missing-cited-path core/docs/broken.md:7: cited path 'core/nope.dart' does not exist",
      "missing-doc-link core/docs/broken.md:7: linked document '../nope.md' does not exist (resolved to core/nope.md)",
      "bad-anchor core/docs/broken.md:7: link 'viewers/README.md#no-such': no heading '#no-such' in core/docs/viewers/README.md",
      "citation-escape core/docs/broken.md:7: link '../../../x.md' escapes the repository",
      "bad-anchor core/docs/broken.md:7: link '#missing': no heading '#missing' in core/docs/broken.md",
    ]);
  });
});

describe('visibility ceiling (review 2 #5)', () => {
  it('refuses a wider visibility than the internal folder allows, accepts narrowing elsewhere', () => {
    expect(homeErrors('core/docs/knowledge-graph/internal/people/bob.yaml', 'id: P:bob\nname: Bob\nemail: b@x\nvisibility: public\n'))
      .toEqual(["visibility-ceiling core/docs/knowledge-graph/internal/people/bob.yaml:4: visibility 'public' exceeds the ceiling 'internal' of core/docs/knowledge-graph/internal/; a home there can only be internal"]);
    expect(homeErrors('core/docs/knowledge-graph/internal/people/bob.yaml', 'id: P:bob\nname: Bob\nemail: b@x\nvisibility: internal\n')).toEqual([]);
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}visibility: internal\n---\n# P\n`)).toEqual([]);
  });
});

describe('hierarchy and edge instances (review 2 #7)', () => {
  it('checks feature roots from schema.yaml, concept depth and the area home of a level-3 feature', () => {
    expect(homeErrors('core/docs/broken.md', '---\nfeature: bogus/thing\nowner: askalkin\n---\n# B\n'))
      .toEqual([expect.stringMatching(/^bad-root core\/docs\/broken.md:2: ~bogus\/thing: feature root 'bogus' is not one of access, transform, .*\(schema.yaml feature_roots\)$/)]);
    expect(homeErrors('core/docs/knowledge-graph/concepts/deep.yaml', 'id: C:a/b/c\nname: Deep\n'))
      .toEqual(['too-deep core/docs/knowledge-graph/concepts/deep.yaml:1: ~C:a/b/c: a concept id has at most two segments', expect.stringMatching(/^stub/)].filter((e) => typeof e === 'string'));
    expect(homeErrors('core/docs/broken.md', '---\nfeature: govern/permissions/roles\n---\n# R\n'))
      .toEqual(['missing-area core/docs/broken.md:2: ~govern/permissions/roles needs its area home ~govern/permissions: a level-2 feature must exist and carry the owner; none found']);
  });

  it('enforces cardinality one, same_type, acyclic and abstract-edge keys on authored edges', () => {
    const repo = makeRepo();
    const kgRoot = path.join(repo, KG_DIR);
    fs.appendFileSync(path.join(kgRoot, 'edges/tracked-in.yaml'), 'cardinality: one\n');
    write(repo, 'core/docs/broken.md', `${PERMISSIONS}tickets: [GROK-1, GROK-2]\n---\n# P\n`);
    fs.writeFileSync(path.join(kgRoot, 'edges/mentors.yaml'), 'type: mentors\nfrom: Person\nto: Person\nkey: mentors\nsame_type: true\nderived_by: [annotation]\ndescription: Who mentors whom.\n');
    write(repo, 'core/docs/knowledge-graph/internal/people/bob.yaml', 'id: P:bob\nname: Bob\nemail: b@x\nmentors: [askalkin]\n');
    write(repo, 'public/help/visualize/viewers/scatter-plot.md', fs.readFileSync(path.join(repo, 'public/help/visualize/viewers/scatter-plot.md'), 'utf8')
      .replace('status: active\n', 'status: active\nsuperseded_by: [visualize/viewers/old-scatter]\n'));
    fs.appendFileSync(path.join(kgRoot, 'edges/evidences.yaml'), 'key: evidence\n');
    write(repo, 'core/docs/evidence.md', `${PERMISSIONS}evidence: [visualize]\n---\n# E\n`.replace('govern/permissions', 'govern/evidence'));
    const {system, homes} = load(repo);
    expect(system.errors).toEqual([]);
    expect(homes.errors.map((e) => `${e.code} ${e.file}: ${e.message}`).sort()).toEqual([
      "acyclic core/docs/viewers/old-scatter.md: supersedes is acyclic but forms a cycle: ~visualize/viewers/scatter-plot -> ~visualize/viewers/old-scatter -> ~visualize/viewers/scatter-plot",
      'bad-edge core/docs/evidence.md: evidence: spells the abstract edge evidences; abstract edges cannot be authored',
      'cardinality core/docs/broken.md: tickets: tracked-in has cardinality one, 2 targets given',
      'cardinality core/docs/viewers/README.md: tickets: tracked-in has cardinality one, 2 targets given',
      'same-type core/docs/knowledge-graph/internal/people/bob.yaml: mentors[0]: mentors is same_type; ~P:askalkin is a developer, this home is a person',
    ].sort());
  });

  it('resolves a reference with an #anchor against the target home\'s headings', () => {
    const tombstone = (target: string) => `---\nfeature: visualize/viewers/older\nstatus: removed\nsuperseded_by: [${target}]\n---\n# Older\n`;
    expect(homeErrors('core/docs/viewers/older.md', tombstone('visualize/viewers/scatter-plot#usage'))).toEqual([]);
    expect(homeErrors('core/docs/viewers/older.md', tombstone('visualize/viewers/scatter-plot#nope')))
      .toEqual(["unresolved-ref core/docs/viewers/older.md:4: superseded_by[0]: 'visualize/viewers/scatter-plot#nope': no heading '#nope' in public/help/visualize/viewers/scatter-plot.md"]);
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}concepts: [dataframe#usage]\n---\n# P\n`))
      .toEqual(["unresolved-ref core/docs/broken.md:4: concepts[0]: 'dataframe#usage': ~C:dataframe is a YAML record and has no headings, so '#usage' cannot resolve"]);
  });
});

describe('duplicate edge-key items', () => {
  it('warns when an edge-key list names the same target twice', () => {
    const repo = makeRepo();
    write(repo, 'core/docs/broken.md', `${PERMISSIONS}concepts: [dataframe, {to: dataframe, role: central}]\n---\n# P\n`);
    const {homes} = load(repo);
    expect(homes.errors).toEqual([]);
    expect(homes.warnings.map((w) => `${w.code} ${w.file}:${w.line}: ${w.message}`))
      .toContain("duplicate-item core/docs/broken.md:4: concepts[1]: 'dataframe' is listed twice under concepts:");
  });
});

describe('member parser shapes and defaults (review 2 #8)', () => {
  it('rejects unions of lists and mixed unions, accepts the parenthesized form', () => {
    expect(parseMember('m?', 'Person[] | Team[]', {types: TYPES}).error).toBe("m: a union of lists is not allowed; write '(Person | Team)[]' instead of 'Person[] | Team[]'");
    expect(parseMember('m?', 'string | string[]').error).toMatch(/may not mix lists and single values/);
    expect(parseMember('m?', '(Person | Team)[]', {types: TYPES}).member).toMatchObject({kind: 'ref', refs: ['person', 'team'], list: true});
  });

  it('validates scalar defaults with the scalar rules, not only typeof', () => {
    expect(parseMember('v?', "Version = '1.x'").error).toMatch(/default "1.x": expected a version like 1.28.0/);
    expect(parseMember('d?', "Date = '2026-13-01'").error).toMatch(/default "2026-13-01": "2026-13-01" is not a calendar date/);
    expect(parseMember('u?', "Url = 'nope'").error).toMatch(/expected an absolute URL/);
    expect(parseMember('p?', "Provenance = 'guess'").error).toMatch(/is not one of 'annotation'/);
    expect(parseMember('d?', "Date = '2026-02-28'").member).toMatchObject({default: '2026-02-28'});
  });

  it('rejects x and x? declared together, and Member outside schema.yaml', () => {
    expect(typeErrors((kg) => append(kg, 'nodes/team.yaml', '  lead:          Person\n')))
      .toEqual([expect.stringMatching(/team\.lead: declared twice \(as 'lead' and 'lead\?'\)/)]);
    expect(parseMember('m?', 'Member').error).toBe("m: 'Member' is the schema's own vocabulary and not allowed in a type file");
    expect(parseMember('m?', 'Member', {schema: true}).member).toMatchObject({scalar: 'Member'});
  });
});

describe('type-file namespaces (review 2 #9)', () => {
  it('reserves edge row fields, rejects keys that shadow node properties, requires prefixes and known inherit names', () => {
    expect(typeErrors((kg) => append(kg, 'edges/covers.yaml', '  from?:         string\n')))
      .toEqual(['edge covers.from: reserved field name, part of every edge row (from, to, type, id)']);
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'edges/owns.yaml'), 'type: owns\nfrom: Feature\nto: Person\nkey: owner\nderived_by: [annotation]\ndescription: Shadows the base owner.\n')))
      .toEqual(["edge key 'owner' collides with the node property 'owner' of node; a frontmatter key can only mean one of them"]);
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'nodes/partner.yaml'), 'type: partner\nextends: actor\nauthored: true\ndescription: No prefix anywhere.\n')))
      .toEqual(['node partner: authored but no type in its chain declares a prefix, so its ids cannot be written']);
    expect(typeErrors((kg) => append(kg, 'nodes/team.yaml', 'inherit: [lead, bogus]\n')))
      .toEqual(["node team: inherit names 'bogus', which is not a member"]);
    expect(typeErrors((kg) => append(kg, 'edges/part-of.yaml', 'inherits: [status, nope]\n')))
      .toEqual(['Feature', 'Concept', 'Scenario'].map((t) => `edge part-of: inherits names 'nope', which is not a member of ${t}`));
    expect(typeErrors((kg) => fs.writeFileSync(path.join(kg, 'nodes/doc--page.yaml'), 'type: doc--page\nextends: artifact\ndescription: Collides.\n')))
      .toEqual(["node types 'doc--page' and 'doc-page' both become DocPage in TypeScript"]);
  });
});

describe('kg.d.ts as an API (review 2 #14)', () => {
  it('escapes literals, carries the build fields and compiles under strict TypeScript', () => {
    const repo = makeRepo();
    fs.appendFileSync(path.join(repo, KG_DIR, 'nodes/team.yaml'), "  mood?:         \"'don\\\\'t' | 'ok'\"\n");
    const system = loadTypeSystem(path.join(repo, KG_DIR));
    expect(system.errors).toEqual([]);
    const dts = generateDts(system);
    expect(dts).toContain("mood?: 'don\\'t' | 'ok';");
    const project = new Project({useInMemoryFileSystem: true, compilerOptions: {strict: true, noEmit: true}});
    project.createSourceFile('kg.d.ts', dts);
    expect(project.getPreEmitDiagnostics().map((d) => d.getMessageText())).toEqual([]);
  });
});

describe('report boundaries (review 2 #15)', () => {
  it('turns unreadable files, unclosed frontmatter and malformed YAML into coded errors', () => {
    const repo = makeRepo();
    write(repo, 'core/docs/open.md', '---\nfeature: govern/x\n# never closed\n');
    write(repo, 'core/docs/bad.md', '---\ntitle: [unclosed\n---\nText.\n');
    expect(load(repo).homes.errors.map((e) => `${e.code} ${e.file}:${e.line}: ${e.message}`)).toEqual([
      expect.stringMatching(/^yaml-error core\/docs\/bad\.md:\d+: YAML error: /),
      'unclosed-frontmatter core/docs/open.md:1: frontmatter is not closed by a --- line',
    ]);
    const system = loadTypeSystem(path.join(repo, KG_DIR));
    expect(loadHomes(system, repo, ['core/docs/vanished.md']).errors)
      .toEqual([{file: 'core/docs/vanished.md', code: 'unreadable', message: expect.stringMatching(/^cannot read: ENOENT/)}]);
  });

  it('applies the discovery excludes when expanding a code: glob', () => {
    const repo = makeRepo();
    write(repo, 'public/packages/Chem/node_modules/dep/index.ts', 'export {};\n');
    expect(homeErrors('core/docs/broken.md', `${PERMISSIONS}code: [public/packages/Chem/**/*.ts]\n---\n# P\n`))
      .toEqual(["missing-path core/docs/broken.md:4: code[0]: no file matches 'public/packages/Chem/**/*.ts'"]);
    void repo;
  });
});
