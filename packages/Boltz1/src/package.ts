/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BoltzService} from './utils/boltz-service';
import {getInitializedService} from './utils/boltz-api-service';
import {ADME_MODEL, BOLTZ_API_CONFIG_PATH, STRUCTURE_AND_BINDING_MODEL} from './utils/boltz-api-constants';
import {scatterRows, toDataFrame} from './utils/boltz-api-utils';
import {Boltz1AppView, openBoltzDemo} from './demo/boltz-app';

export * from './package.g';
export const _package = new DG.Package();

async function getBoltzApiConfigFolders(taskType: string): Promise<string[]> {
  const files = await grok.dapi.files.list(`${BOLTZ_API_CONFIG_PATH}/${taskType}`, true);
  return files.filter((f) => f.isDirectory).map((f) => f.name);
}

async function readBoltzApiConfig(taskType: string, config: string): Promise<any> {
  return JSON.parse(await grok.dapi.files.readAsText(`${BOLTZ_API_CONFIG_PATH}/${taskType}/${config}/config.json`));
}

export class PackageFunctions {
  @grok.decorators.func()
  static info(): void {
    grok.shell.info(_package.webRoot);
  }

  @grok.decorators.func()
  static async getBoltzConfigFolders(): Promise<string[]> {
    return await BoltzService.getBoltzConfigFolders();
  }

  @grok.decorators.func({
    meta: {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *',
    },
  })
  static async runBoltz(
    config: string,
    msa: string,
  ): Promise<string> {
    return await BoltzService.runBoltz(config, msa);
  }

  @grok.decorators.func({
    'name': 'Boltz',
    'top-menu': 'Bio | Folding | Boltz...',
    'meta': {vectorFunc: 'true'},
    'outputs': [{name: 'result', type: 'dataframe', options: {action: 'join(table)'}}],
  })
  static async folding(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Macromolecule'}}) sequences: DG.Column,
  ): Promise<DG.DataFrame> {
    return await BoltzService.folding(table, sequences);
  }

  @grok.decorators.func({
    'name': 'Boltz',
    'top-menu': 'Chem | Docking | Boltz...',
    'meta': {vectorFunc: 'true'},
    'outputs': [{name: 'result', type: 'dataframe', options: {action: 'join(table)'}}],
  })
  static async docking(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) ligands: DG.Column,
    @grok.decorators.param({options: {choices: 'Boltz1:getBoltzConfigFolders', description: '\'Folder with config files for docking\''}}) config: string,
  ): Promise<DG.DataFrame> {
    return await BoltzService.docking(table, ligands, config);
  }

  @grok.decorators.panel({
    name: 'Boltz-1',
    outputs: [
      {name: 'result', type: 'widget'},
    ],
    condition: 'Boltz1:isApplicableBoltz(molecule)',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static async boltzWidget(
    @grok.decorators.param({options: {semType: 'Molecule3D'}}) molecule: DG.SemanticValue,
  ): Promise<DG.Widget<any> | null> {
    return await BoltzService.boltzWidget(molecule);
  }

  @grok.decorators.func()
  static isApplicableBoltz(
    molecule: string,
  ): boolean {
    return molecule.includes('binding_confidence');
  }

  @grok.decorators.app({
    name: 'Boltz-1',
    browsePath: 'Bio',
  })
  static async boltz1App(): Promise<DG.ViewBase> {
    return new Boltz1AppView().getView();
  }

  @grok.decorators.func({
    name: 'Boltz Demo',
    description: 'WDR5 protein-ligand complexes with Boltz-predicted affinity',
    meta: {demoPath: 'Bioinformatics | Boltz', isDemoDashboard: 'true'},
  })
  static async demoBoltz(): Promise<void> {
    await openBoltzDemo();
  }

  // --- Hosted API: config-listing helpers (used as choices providers) ---

  @grok.decorators.func()
  static async getBoltzStructureBindingConfigs(): Promise<string[]> {
    return getBoltzApiConfigFolders('structure-binding');
  }

  @grok.decorators.func()
  static async getBoltzSmDesignConfigs(): Promise<string[]> {
    return getBoltzApiConfigFolders('sm-design');
  }

  @grok.decorators.func()
  static async getBoltzSmScreenConfigs(): Promise<string[]> {
    return getBoltzApiConfigFolders('sm-screen');
  }

  @grok.decorators.func()
  static async getBoltzProteinDesignConfigs(): Promise<string[]> {
    return getBoltzApiConfigFolders('protein-design');
  }

  @grok.decorators.func()
  static async getBoltzProteinScreenConfigs(): Promise<string[]> {
    return getBoltzApiConfigFolders('protein-screen');
  }

  // --- Hosted API: main functions ---

  @grok.decorators.func({
    'meta': {vectorFunc: 'true'},
    'outputs': [{name: 'result', type: 'dataframe', options: {action: 'join(table)'}}],
  })
  static async boltzStructureAndBinding(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) ligands: DG.Column,
    @grok.decorators.param({options: {choices: 'Boltz1:getBoltzStructureBindingConfigs'}}) config: string,
  ): Promise<DG.DataFrame> {
    const svc = await getInitializedService();
    const cfg = await readBoltzApiConfig('structure-binding', config);
    const ligandChain = cfg.ligand_chain_id ?? 'B';
    const jobs = await Promise.all((ligands.toList() as string[]).map((smiles) => svc.predictStructureAndBinding({
      model: STRUCTURE_AND_BINDING_MODEL,
      input: {
        entities: [
          {type: 'protein', value: cfg.protein.value, chain_ids: cfg.protein.chain_ids},
          {type: 'ligand_smiles', value: smiles, chain_ids: [ligandChain]},
        ],
        binding: {type: 'ligand_protein_binding', binder_chain_id: ligandChain},
      },
    })));
    const complexes = await Promise.all(jobs.map(async (j) => {
      const url = j.output?.best_sample?.structure?.url;
      return url ? await (await grok.dapi.fetchProxy(url)).text() : '';
    }));
    const result = toDataFrame(jobs.map((j) => ({
      ...j.output?.best_sample?.metrics,
      binding_confidence: j.output?.binding_metrics?.binding_confidence,
      optimization_score: (j.output?.binding_metrics as any)?.optimization_score,
    })));
    const complexCol = DG.Column.fromStrings('Complex', complexes);
    complexCol.semType = DG.SEMTYPE.MOLECULE3D;
    complexCol.setTag(DG.TAGS.CELL_RENDERER, DG.SEMTYPE.MOLECULE3D);
    result.columns.add(complexCol);
    return result;
  }

  @grok.decorators.func({
    'meta': {vectorFunc: 'true'},
    'outputs': [{name: 'result', type: 'dataframe', options: {action: 'join(table)'}}],
  })
  static async boltzAdme(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
  ): Promise<DG.DataFrame> {
    const svc = await getInitializedService();
    const smiles = molecules.toList() as string[];
    const batches = [];
    for (let i = 0; i < smiles.length; i += 128)
      batches.push(smiles.slice(i, i + 128).map((s, j) => ({smiles: s, id: String(i + j)})));
    const jobs = await Promise.all(batches.map((mols) => svc.predictAdme({model: ADME_MODEL, input: {molecules: mols}})));
    const rows = Array.from({length: smiles.length},
      () => ({lipophilicity: NaN, permeability: NaN, solubility: null as string | null}));
    for (const job of jobs)
      for (const mol of (job.output?.molecules ?? []))
        if (mol.status === 'succeeded')
          rows[parseInt(mol.external_id!)] = {
            lipophilicity: mol.adme?.lipophilicity ?? NaN,
            permeability: mol.adme?.permeability ?? NaN,
            solubility: mol.adme?.solubility ?? null,
          };
    return toDataFrame(rows);
  }

  @grok.decorators.func()
  static async boltzDesignSmallMolecules(
    @grok.decorators.param({options: {choices: 'Boltz1:getBoltzSmDesignConfigs'}}) config: string,
    numMolecules: number,
  ): Promise<DG.DataFrame> {
    const svc = await getInitializedService();
    const cfg = await readBoltzApiConfig('sm-design', config);
    const results = await svc.designSmallMolecules({target: cfg.target, num_molecules: numMolecules});
    return toDataFrame(results.map((r) => ({smiles: r.smiles, ...r.metrics})));
  }

  @grok.decorators.func({
    'meta': {vectorFunc: 'true'},
    'outputs': [{name: 'result', type: 'dataframe', options: {action: 'join(table)'}}],
  })
  static async boltzScreenSmallMolecules(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({options: {choices: 'Boltz1:getBoltzSmScreenConfigs'}}) config: string,
  ): Promise<DG.DataFrame> {
    const svc = await getInitializedService();
    const cfg = await readBoltzApiConfig('sm-screen', config);
    const results = await svc.screenSmallMoleculeLibrary({
      target: cfg.target,
      molecules: (molecules.toList() as string[]).map((s, i) => ({smiles: s, id: String(i)})),
    });
    return toDataFrame(scatterRows(table.rowCount, results.map((r) => ({
      index: parseInt(r.external_id ?? '-1'),
      row: {smiles: r.smiles, ...r.metrics},
    }))));
  }

  @grok.decorators.func()
  static async boltzDesignProteins(
    @grok.decorators.param({options: {choices: 'Boltz1:getBoltzProteinDesignConfigs'}}) config: string,
    numProteins: number,
  ): Promise<DG.DataFrame> {
    const svc = await getInitializedService();
    const cfg = await readBoltzApiConfig('protein-design', config);
    const binderChain = cfg.binder_specification?.entities?.[0]?.chain_ids?.[0] ?? 'B';
    const results = await svc.designProteins({
      target: cfg.target,
      binder_specification: cfg.binder_specification,
      num_proteins: numProteins,
    });
    return toDataFrame(results.map((r) => ({
      sequence: r.entities.find((e: any) => e.type === 'protein' && e.chain_ids.includes(binderChain))?.value,
      ...r.metrics,
    })));
  }

  @grok.decorators.func({
    'meta': {vectorFunc: 'true'},
    'outputs': [{name: 'result', type: 'dataframe', options: {action: 'join(table)'}}],
  })
  static async boltzScreenProteins(
    table: DG.DataFrame,
    @grok.decorators.param({options: {semType: 'Macromolecule'}}) proteins: DG.Column,
    @grok.decorators.param({options: {choices: 'Boltz1:getBoltzProteinScreenConfigs'}}) config: string,
  ): Promise<DG.DataFrame> {
    const svc = await getInitializedService();
    const cfg = await readBoltzApiConfig('protein-screen', config);
    const results = await svc.screenProteinLibrary({
      target: cfg.target,
      proteins: (proteins.toList() as string[]).map((seq, i) => ({
        id: String(i),
        entities: [{type: 'protein', value: seq, chain_ids: [cfg.binder_chain_id ?? 'B']}],
      })),
    });
    return toDataFrame(scatterRows(table.rowCount, results.map((r) => ({
      index: parseInt(r.external_id ?? '-1'),
      row: {...r.metrics},
    }))));
  }
}
