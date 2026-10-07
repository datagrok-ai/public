/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HmmerPool} from './pool';
import {allowedChains, numberColumn, numberingFrame, SCHEMES, speciesOption} from './numbering';
import {resolveModels, searchColumn, showProfileSearchDialog} from './search';
import {findDomains, showFindDomainsDialog} from './domains';
import {anarciWidget, domainsWidget} from './widgets';
import type {Scheme} from './hmmer/anarci/types.ts';

export * from './package.g';
export const _package = new DG.Package();

export function pool(): HmmerPool {
  return HmmerPool.get(_package);
}

export class PackageFunctions {
  @grok.decorators.func({
    name: 'anarciNumbering',
    friendlyName: 'ANARCI (HMMER)',
    description: 'ANARCI numbering (IMGT, Kabat, Chothia, Martin, AHo, Wolfguy) of antibody and T-cell receptor chains, computed in the browser by an exact WebAssembly port of HMMER 3.4. Also returns chain, species, V/J germlines and E-values.',
    meta: {role: 'antibodyNumbering'},
  })
  static async anarciNumbering(
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    @grok.decorators.param({type: 'dataframe'}) df: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule'}}) seqCol: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {choices: ['imgt', 'kabat', 'chothia', 'martin', 'aho', 'wolfguy'], initialValue: 'imgt'}}) scheme: string,
  ): Promise<DG.DataFrame> {
    const s = (SCHEMES.includes(scheme as Scheme) ? scheme : 'imgt') as Scheme;
    const {results} = await numberColumn(pool(), seqCol,
      {scheme: s, allow: allowedChains(s), assignGermline: true});
    return numberingFrame(s, results);
  }

  @grok.decorators.func({
    name: 'antibodyGermlines',
    friendlyName: 'Antibody Germlines (ANARCI)',
    description: 'Adds chain type, species, closest V and J germline genes with identities, and the HMMER E-value of antibody and T-cell receptor chains (ANARCI method, computed in the browser).',
    topMenu: 'Bio | Annotate | Germlines and Species (ANARCI)...',
    outputs: [],
  })
  static async antibodyGermlines(
    @grok.decorators.param({options: {caption: 'Table', nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column',
      options: {semType: 'Macromolecule', caption: 'Sequence', nullable: false}}) sequence: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {caption: 'Species', choices: ['human, mouse', 'all', 'human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig', 'alpaca', 'cow'], initialValue: 'human, mouse', description: 'Preferred species (ANARCI default: human and mouse)'}})
      species: string = 'human, mouse',
  ): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('ANARCI germlines...');
    try {
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
    description: 'Searches a sequence column with a profile HMM (Pfam family, Pfam accession or HMM file) using hmmsearch in the browser; adds score, E-value and domain columns and annotates the domains.',
    topMenu: 'Bio | Search | Profile HMM Search...',
  })
  static async profileHmmSearchDialog(): Promise<void> {
    await showProfileSearchDialog(pool(), _package);
  }

  @grok.decorators.func({
    name: 'searchWithHmm',
    description: 'hmmsearch of a sequence column with one profile HMM. Model: HMMER text, Pfam accession(s) (PF07686) or builtin:<index>. Returns score, E-value, domains, significant and regions columns.',
    outputs: [{name: 'result', type: 'dataframe'}],
  })
  static async searchWithHmm(
    @grok.decorators.param({options: {caption: 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule', caption: 'Sequence'}}) sequence: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {caption: 'Model'}}) model: string,
    @grok.decorators.param({type: 'double', options: {caption: 'E-value', initialValue: '10'}}) evalue: number = 10,
    @grok.decorators.param({type: 'bool', options: {caption: 'Gathering cutoff', initialValue: 'false'}}) gathering: boolean = false,
  ): Promise<DG.DataFrame> {
    const models = await resolveModels(pool(), _package, model);
    const index = model.startsWith('builtin:') ? Number(model.slice(8)) : 0;
    const result = await searchColumn(pool(), models, index, sequence, gathering ? {cutoff: 'ga'} : {evalue});
    return DG.DataFrame.fromColumns(result.columns);
  }

  @grok.decorators.func({
    name: 'findDomainsDialog',
    friendlyName: 'Find Domains (HMMER)',
    description: 'Annotates protein domains (built-in Pfam biologics library, Pfam accessions or an HMM file) on a sequence column using hmmscan in the browser.',
    topMenu: 'Bio | Annotate | Find Domains (HMMER)...',
  })
  static async findDomainsDialog(): Promise<void> {
    await showFindDomainsDialog(pool(), _package);
  }

  @grok.decorators.func({
    name: 'findDomains',
    description: 'hmmscan of a sequence column against a model library: builtin, Pfam accession(s) or HMMER text. Returns one row per domain (row, model, accession, from, to, score, evalue).',
    outputs: [{name: 'result', type: 'dataframe'}],
  })
  static async findDomains(
    @grok.decorators.param({options: {caption: 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {semType: 'Macromolecule', caption: 'Sequence'}}) sequence: DG.Column<string>,
    @grok.decorators.param({type: 'string', options: {caption: 'Library', initialValue: 'builtin'}}) library: string = 'builtin',
    @grok.decorators.param({type: 'double', options: {caption: 'E-value', initialValue: '0.01'}}) evalue: number = 0.01,
    @grok.decorators.param({type: 'bool', options: {caption: 'Gathering cutoffs', initialValue: 'true'}}) gathering: boolean = true,
  ): Promise<DG.DataFrame> {
    const models = await resolveModels(pool(), _package, library);
    const {found} = await findDomains(pool(), models, sequence, gathering ? {cutoff: 'ga'} : {evalue, domEvalue: evalue});
    return DG.DataFrame.fromObjects(found) ?? DG.DataFrame.create(0);
  }

  @grok.decorators.panel({
    name: 'Bioinformatics | ANARCI',
    description: 'Chain, species, germlines and IMGT CDRs of an antibody or TCR sequence (ANARCI on HMMER)',
    tags: ['bio', 'widgets', 'panel'],
    meta: {role: 'widgets', domain: 'bio'},
  })
  static anarciPanel(
    @grok.decorators.param({options: {semType: 'Macromolecule', units: 'fasta'}}) sequence: DG.SemanticValue): DG.Widget {
    return anarciWidget(pool(), sequence);
  }

  @grok.decorators.panel({
    name: 'Bioinformatics | Pfam Domains',
    description: 'Pfam domains of a protein sequence (built-in library, gathering cutoffs; HMMER hmmscan)',
    tags: ['bio', 'widgets', 'panel'],
    meta: {role: 'widgets', domain: 'bio'},
  })
  static pfamDomainsPanel(
    @grok.decorators.param({options: {semType: 'Macromolecule', units: 'fasta'}}) sequence: DG.SemanticValue): DG.Widget {
    return domainsWidget(pool(), _package, sequence);
  }

  @grok.decorators.func({
    name: 'hmmerEngineInfo',
    description: 'HMMER engine status: worker count and loaded model databases.',
    outputs: [{name: 'result', type: 'object'}],
  })
  static hmmerEngineInfo(): object {
    return {maxWorkers: HmmerPool.maxWorkers, webRoot: _package.webRoot};
  }
}
