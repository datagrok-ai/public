/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HmmerPool} from './pool';
import type {Scheme} from './hmmer/anarci/types.ts';

// Feature modules load on first use, keeping package.js small (web/budgets.json).
const numbering = () => import('./numbering');
const search = () => import('./search');
const domains = () => import('./domains');

/** A context panel filled asynchronously by a lazily loaded module. */
function lazyWidget(fill: () => Promise<DG.Widget>): DG.Widget {
  const root = ui.div([ui.loader()]);
  fill().then((w) => root.replaceChildren(w.root)).catch((e) =>
    root.replaceChildren(ui.divText(`${e instanceof Error ? e.message : e}`, {style: {color: 'var(--red-3)'}})));
  return new DG.Widget(root);
}

export * from './package.g';
export const _package = new DG.Package();

export function pool(): HmmerPool {
  return HmmerPool.get(_package);
}

export class PackageFunctions {
  @grok.decorators.func({
    name: 'anarciNumbering',
    friendlyName: 'ANARCI (HMMER)',
    description: 'ANARCI numbering (IMGT, Kabat, Chothia, Martin, AHo, Wolfguy) of antibody and T-cell receptor chains, computed in the browser. Also returns chain, species, V/J germlines and E-values. ANARCI and HMMER are reimplemented here (TypeScript; Rust compiled to WebAssembly), not the original programs; results are identical.',
    meta: {role: 'antibodyNumbering'},
  })
  static async anarciNumbering(
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    @grok.decorators.param({type: 'dataframe'}) df: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule'}}) seqCol: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {choices: ['imgt', 'kabat', 'chothia', 'martin', 'aho', 'wolfguy'], initialValue: 'imgt'}}) scheme: string,
    @grok.decorators.param({type: 'string', options: {optional: true, caption: 'Species', choices: ['human and mouse', 'all', 'human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig', 'alpaca', 'cow'], initialValue: 'human and mouse', description: 'Preferred species of the germline models (ANARCI default: human and mouse, other species only when they score clearly better)'}})
      species: string = 'human and mouse',
    @grok.decorators.param({type: 'double', options: {optional: true, caption: 'Bit score threshold', initialValue: '80', description: 'Minimum HMMER bit score of a domain (ANARCI default 80)'}})
      bitScoreThreshold: number = 80,
    @grok.decorators.param({type: 'string', options: {optional: true, caption: 'Domain', choices: ['first', 'heavy', 'light'], initialValue: 'first', description: 'Which domain of multi-domain chains (scFv, diabodies) to number: the N-terminal one, the first heavy (TCR beta/delta) or the first light (TCR alpha/gamma) domain'}})
      domain: string = 'first',
  ): Promise<DG.DataFrame> {
    const {allowedChains, domainChoice, numberColumn, numberingFrame, SCHEMES, speciesOption} = await numbering();
    const s = (SCHEMES.includes(scheme as Scheme) ? scheme : 'imgt') as Scheme;
    const threshold = Number.isFinite(bitScoreThreshold) ? bitScoreThreshold : 80;
    const {results} = await numberColumn(pool(), seqCol, {scheme: s, allow: allowedChains(s), assignGermline: true,
      allowedSpecies: speciesOption(species), bitScoreThreshold: threshold});
    return numberingFrame(s, results, domainChoice(domain));
  }

  @grok.decorators.func({
    name: 'antibodyGermlines',
    friendlyName: 'Antibody Germlines (ANARCI)',
    description: 'Adds chain type, species, closest V and J germline genes with identities, and the HMMER E-value of antibody and T-cell receptor chains (ANARCI method, computed in the browser). ANARCI and HMMER are reimplemented here (TypeScript; Rust compiled to WebAssembly), not the original programs; results are identical.',
    topMenu: 'Bio | Annotate | Germlines and Species (ANARCI)...',
    outputs: [],
  })
  static async antibodyGermlines(
    @grok.decorators.param({options: {caption: 'Table', nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column',
      options: {semType: 'Macromolecule', caption: 'Sequence', nullable: false}}) sequence: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {caption: 'Species', choices: ['human and mouse', 'all', 'human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig', 'alpaca', 'cow'], initialValue: 'human and mouse', description: 'Preferred species (ANARCI default: human and mouse)'}})
      species: string = 'human and mouse',
  ): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('ANARCI germlines...');
    try {
      const {numberColumn, numberingFrame, speciesOption} = await numbering();
      const {results} = await numberColumn(pool(), sequence,
        {scheme: 'imgt', assignGermline: true, allowedSpecies: speciesOption(species)});
      const frame = numberingFrame('imgt', results);
      const prefix = `${sequence.name} `;
      for (const name of ['chain', 'species', 'v_gene', 'v_identity', 'j_gene', 'j_identity', 'evalue', 'bitscore']) {
        const col = frame.getCol(name);
        col.name = table.columns.getUnusedName(prefix + name.replace('_', ' '));
        table.columns.add(col);
      }
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func({
    name: 'profileHmmSearchDialog',
    friendlyName: 'Profile HMM Search',
    description: 'Searches a sequence column with a profile HMM (Pfam family, Pfam accession or HMM file) using hmmsearch in the browser; adds score, E-value and domain columns and annotates the domains. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    topMenu: 'Bio | Search | Profile HMM Search...',
  })
  static async profileHmmSearchDialog(): Promise<void> {
    await (await search()).showProfileSearchDialog(pool(), _package);
  }

  @grok.decorators.func({
    name: 'searchWithHmm',
    description: 'hmmsearch of a sequence column with one profile HMM. Model: HMMER text, Pfam accession(s) (PF07686) or builtin:<index>. Returns score, E-value, domains, significant and regions columns. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    outputs: [{name: 'result', type: 'dataframe'}],
  })
  static async searchWithHmm(
    @grok.decorators.param({options: {caption: 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule', caption: 'Sequence'}}) sequence: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {caption: 'Model'}}) model: string,
    @grok.decorators.param({type: 'double', options: {caption: 'E-value', initialValue: '10'}}) evalue: number = 10,
    @grok.decorators.param({type: 'bool', options: {caption: 'Gathering cutoff', initialValue: 'false'}}) gathering: boolean = false,
  ): Promise<DG.DataFrame> {
    const {resolveModels, searchColumn} = await search();
    const models = await resolveModels(pool(), _package, model);
    const index = model.startsWith('builtin:') ? Number(model.slice(8)) : 0;
    const result = await searchColumn(pool(), models, index, sequence, gathering ? {cutoff: 'ga'} : {evalue});
    return DG.DataFrame.fromColumns(result.columns);
  }

  @grok.decorators.func({
    name: 'alignToProfile',
    friendlyName: 'Profile HMM (hmmalign)',
    description: 'Aligns sequences to one profile HMM with hmmalign (HMMER 3.4, in the browser): the best-matching Pfam family of the built-in library, a Pfam family or accession, the best antibody/TCR germline model (ANARCI) or an HMM file. Columns are model positions and insertions between them. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    meta: {role: 'sequenceMSA', acceptsCanonical: 'true'},
    outputs: [{name: 'result', type: 'column'}],
  })
  static async alignToProfile(
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule'}}) sequences: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {optional: true, caption: 'Profile', choices: ['Best Pfam family (auto)', 'Pfam family or accession', 'Antibody/TCR germline (auto)', 'HMM file'], initialValue: 'Best Pfam family (auto)', description: 'The profile HMM to align to'}})
      profile: string = 'Best Pfam family (auto)',
    @grok.decorators.param({type: 'string', options: {optional: true, nullable: true, caption: 'Family', description: 'Pfam family of the built-in library (V-set) or any Pfam accession (PF07686), for Profile = Pfam family or accession'}})
      family: string = '',
    @grok.decorators.param({type: 'file', options: {optional: true, nullable: true, caption: 'HMM file', description: 'A file with one profile HMM (HMMER2/3 text or .h3m), for Profile = HMM file'}})
      hmmFile: DG.FileInfo | null = null,
    @grok.decorators.param({type: 'bool', options: {optional: true, caption: 'Trim', initialValue: 'false', description: 'Trim residues outside the model\'s alignment (hmmalign --trim)'}})
      trim: boolean = false,
  ): Promise<DG.Column<string>> {
    const [{alignColumn, chooseProfile}, {columnResidues}] = await Promise.all([import('./align'), import('./sequences')]);
    const named = (await columnResidues(sequences)).residues.map((residues, i) => ({name: String(i), residues}));
    const {models, index} = await chooseProfile(pool(), _package, profile || 'Best Pfam family (auto)', family,
      hmmFile, named);
    return alignColumn(pool(), sequences, models, index, trim === true, named);
  }

  @grok.decorators.func({
    name: 'hmmWebLogoProfile',
    friendlyName: 'Profile HMM Logo Heights',
    description: 'Profiles for Bio\'s WebLogo (its Profile property): hmmlogo heights in bits of a Pfam family of the built-in library (V-set), a Pfam accession (PF07686), an ANARCI germline model (germline:human_H) or a profile built in the table. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    meta: {role: 'webLogoProfile'},
    outputs: [{name: 'result', type: 'object'}],
  })
  static async hmmWebLogoProfile(
    @grok.decorators.param({type: 'dataframe', options: {nullable: true}}) table: DG.DataFrame | null,
    @grok.decorators.param({type: 'string', options: {nullable: true}}) column: string,
    @grok.decorators.param({type: 'string'}) profile: string,
  ): Promise<object> {
    const {webLogoProfile} = await import('./logo');
    return webLogoProfile(pool(), _package, table, column ?? '', profile);
  }

  @grok.decorators.func({
    name: 'profileLogoDialog',
    friendlyName: 'Profile HMM Logo',
    description: 'Adds a WebLogo of a profile HMM (hmmlogo heights) to the current table view. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    topMenu: 'Bio | Analyze | Profile HMM Logo...',
  })
  static async profileLogoDialog(): Promise<void> {
    await (await import('./logo')).showLogoDialog(pool(), _package);
  }

  @grok.decorators.func({
    name: 'buildProfileDialog',
    friendlyName: 'Build Profile HMM',
    description: 'Builds a profile HMM (hmmbuild, in the browser) from aligned rows of a column, then searches the table with it, annotates the hits and shows its logo. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    topMenu: 'Bio | Search | Build Profile HMM...',
  })
  static async buildProfileDialog(): Promise<void> {
    await (await import('./build')).showBuildProfileDialog(pool());
  }

  @grok.decorators.func({
    name: 'buildProfile',
    description: 'hmmbuild of the rows of an aligned sequence column (all, or the selected ones). Returns the HMMER3 model text, identical to C hmmbuild -n name; search with it via searchWithHmm. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    outputs: [{name: 'result', type: 'string'}],
  })
  static async buildProfile(
    @grok.decorators.param({options: {caption: 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule', caption: 'Alignment'}}) alignment: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {caption: 'Name'}}) name: string,
    @grok.decorators.param({type: 'bool', options: {optional: true, caption: 'Selected rows only', initialValue: 'false'}}) selectedOnly: boolean = false,
  ): Promise<string> {
    const rows = selectedOnly ? Array.from(table.selection.getSelectedIndexes()) :
      Array.from({length: table.rowCount}, (_, i) => i);
    return (await (await import('./build')).buildProfile(pool(), alignment, rows, name)).text;
  }

  @grok.decorators.func({
    name: 'homologySearchDialog',
    friendlyName: 'Homology Search (phmmer)',
    description: 'Searches a protein sequence column with one query sequence (a row or typed) using phmmer in the browser; adds bit score, E-value and domain columns and annotates the matches. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    topMenu: 'Bio | Search | Homology Search (phmmer)...',
  })
  static async homologySearchDialog(): Promise<void> {
    const {showHomologySearchDialog} = await import('./homology');
    await showHomologySearchDialog(pool());
  }

  @grok.decorators.func({
    name: 'homologySearch',
    description: 'phmmer of one query protein sequence against a sequence column (Z = non-empty rows). Returns score, E-value, domains, significant and regions columns. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    outputs: [{name: 'result', type: 'dataframe'}],
  })
  static async homologySearch(
    @grok.decorators.param({options: {caption: 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule', caption: 'Sequence'}}) sequence: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {caption: 'Query'}}) query: string,
    @grok.decorators.param({type: 'double', options: {optional: true, caption: 'E-value', initialValue: '10'}}) evalue: number = 10,
    @grok.decorators.param({type: 'string', options: {optional: true, caption: 'Matrix', choices: ['BLOSUM62', 'BLOSUM45', 'BLOSUM50', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM120', 'PAM240'], initialValue: 'BLOSUM62'}}) matrix: string = 'BLOSUM62',
  ): Promise<DG.DataFrame> {
    const {homologySearch} = await import('./homology');
    const {result} = await homologySearch(pool(), sequence, query, {evalue: evalue ?? 10, matrix: matrix || 'BLOSUM62'});
    return DG.DataFrame.fromColumns(result.columns);
  }

  @grok.decorators.func({
    name: 'findDomainsDialog',
    friendlyName: 'Find Domains (HMMER)',
    description: 'Annotates protein domains (built-in Pfam biologics library, Pfam accessions or an HMM file) on a sequence column using hmmscan in the browser. Overlapping matches of one Pfam clan are resolved as Pfam does (this package\'s implementation, checked against pfam_scan.pl). HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    topMenu: 'Bio | Annotate | Find Domains (HMMER)...',
  })
  static async findDomainsDialog(): Promise<void> {
    await (await domains()).showFindDomainsDialog(pool(), _package);
  }

  @grok.decorators.func({
    name: 'findDomains',
    description: 'hmmscan of a sequence column against a model library: builtin, Pfam accession(s) or HMMER text. Returns one row per domain (row, model, accession, from, to, score, evalue, clan); overlapping domains of a Pfam clan are resolved as Pfam does. HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    outputs: [{name: 'result', type: 'dataframe'}],
  })
  static async findDomains(
    @grok.decorators.param({options: {caption: 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule', caption: 'Sequence'}}) sequence: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {caption: 'Library', initialValue: 'builtin'}}) library: string = 'builtin',
    @grok.decorators.param({type: 'double', options: {caption: 'E-value', initialValue: '0.01'}}) evalue: number = 0.01,
    @grok.decorators.param({type: 'bool', options: {caption: 'Gathering cutoffs', initialValue: 'true'}}) gathering: boolean = true,
    @grok.decorators.param({type: 'bool', options: {optional: true, caption: 'Clan competition', initialValue: 'true', description: 'Of overlapping domains from the same Pfam clan keep the most significant (Pfam families only)'}}) clanCompetition: boolean = true,
  ): Promise<DG.DataFrame> {
    const [{resolveModels}, {findDomains}] = await Promise.all([search(), domains()]);
    const models = await resolveModels(pool(), _package, library);
    const {found} = await findDomains(pool(), models, sequence,
      gathering ? {cutoff: 'ga'} : {evalue, domEvalue: evalue}, clanCompetition !== false);
    return DG.DataFrame.fromObjects(found) ?? DG.DataFrame.create(0);
  }

  @grok.decorators.panel({
    name: 'Bioinformatics | ANARCI',
    description: 'Chain, species, germlines and IMGT CDRs of an antibody or TCR sequence (ANARCI method). ANARCI and HMMER are reimplemented here (TypeScript; Rust compiled to WebAssembly), not the original programs; results are identical.',
    tags: ['bio', 'widgets', 'panel'],
    meta: {role: 'widgets', domain: 'bio'},
  })
  static anarciPanel(
    @grok.decorators.param({options: {semType: 'Macromolecule', units: 'fasta'}}) sequence: DG.SemanticValue): DG.Widget {
    return lazyWidget(async () => (await import('./widgets')).anarciWidget(pool(), sequence));
  }

  @grok.decorators.panel({
    name: 'Bioinformatics | Pfam Domains',
    description: 'Pfam domains of a protein sequence (built-in library, gathering cutoffs; hmmscan). HMMER 3.4 is reimplemented here in Rust and compiled to WebAssembly (Rusty-HMMER), not the original HMMER; results are identical.',
    tags: ['bio', 'widgets', 'panel'],
    meta: {role: 'widgets', domain: 'bio'},
  })
  static pfamDomainsPanel(
    @grok.decorators.param({options: {semType: 'Macromolecule', units: 'fasta'}}) sequence: DG.SemanticValue): DG.Widget {
    return lazyWidget(async () => (await import('./widgets')).domainsWidget(pool(), _package, sequence));
  }

  @grok.decorators.func({
    name: 'hmmerEngineInfo',
    description: 'Engine status (Rusty-HMMER, a Rust/WebAssembly reimplementation of HMMER 3.4): worker count and loaded model databases.',
    outputs: [{name: 'result', type: 'object'}],
  })
  static hmmerEngineInfo(): object {
    return {maxWorkers: HmmerPool.maxWorkers, webRoot: _package.webRoot};
  }
}
