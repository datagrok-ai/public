/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { BoltzService } from './utils/boltz-service';
import { Boltz1AppView } from './demo/boltz-app';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.func({
    description: 'Show package web root',
  })
  static info(): void {
    grok.shell.info(_package.webRoot);
  }

  @grok.decorators.func()
  static async getBoltzConfigFolders(): Promise<string[]> {
    return await BoltzService.getBoltzConfigFolders();
  }

  @grok.decorators.func({
    meta: {
      cache: 'all',
      'cache.invalidateOn': '0 0 1 * *',
    }
  })
  static async runBoltz(
    @grok.decorators.param({name: 'config', type: DG.TYPE.STRING}) config: string,
    @grok.decorators.param({name: 'msa', type: DG.TYPE.STRING}) msa: string,
  ): Promise<string> {
    return await BoltzService.runBoltz(config, msa);
  }

  @grok.decorators.func({
    name: 'Folding',
    'top-menu': 'Bio | Folding | Boltz-1...'
  })
  static async folding(
    @grok.decorators.param({name: 'df', type: DG.TYPE.DATA_FRAME}) df: DG.DataFrame,
    @grok.decorators.param({name: 'sequences', type: DG.TYPE.COLUMN, options: {semType: 'Macromolecule'}}) sequences: DG.Column,
  ): Promise<DG.DataFrame> {
    return await BoltzService.folding(df, sequences);
  }

  @grok.decorators.func({
    name: 'Docking',
    'top-menu': 'Chem | Docking | Boltz-1...'
  })
  static async docking(
    @grok.decorators.param({name: 'df', type: DG.TYPE.DATA_FRAME}) df: DG.DataFrame,
    @grok.decorators.param({name: 'molecules', type: DG.TYPE.COLUMN, options: {semType: 'Molecule'}}) molecules: DG.Column,
    @grok.decorators.param({name: 'config', type: DG.TYPE.STRING, options: {choices: 'Boltz1: getBoltzConfigFolders', caption: 'Folder with config files for docking'}}) config: string,
  ): Promise<DG.DataFrame> {
    return await BoltzService.docking(df, molecules, config);
  }

  @grok.decorators.panel({
    name: 'Boltz-1',
    tags: ['chem', 'widgets'],
    outputs: [
      {name: 'result', type: 'widget'},
    ],
    condition: 'Boltz1:isApplicableBoltz(molecule)',
  })
  static async boltzWidget(
    @grok.decorators.param({name: 'molecule', type: DG.TYPE.SEMANTIC_VALUE, options: {semType: 'Molecule3D'}}) molecule: DG.SemanticValue,
  ): Promise<DG.Widget<any> | null> {
    return await BoltzService.boltzWidget(molecule);
  }

  @grok.decorators.func()
  static isApplicableBoltz(
    @grok.decorators.param({name: 'molecule', type: DG.TYPE.STRING}) molecule: string,
  ): boolean {
    return molecule.includes('confidence_score');
  }

  @grok.decorators.app({
    name: 'Boltz-1',
    browsePath: 'Bio',
  })
  static async boltz1App(): Promise<DG.ViewBase> {
    return new Boltz1AppView().getView();
  }
}