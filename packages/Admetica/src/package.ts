/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  convertLD50,
  getModelsSingle,
  performChemicalPropertyPredictions,
  setProperties,
} from './utils/admetica-utils';
import { properties } from './utils/admetica-utils';
import { AdmeticaBaseEditor } from './utils/admetica-editor';
import { Model, Subgroup } from './utils/constants';
import { AdmeticaViewApp } from './utils/admetica-app';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  // @grok.decorators.init()
  // static async init() { }

  @grok.decorators.func()
  static info() {
    grok.shell.info(_package.webRoot);
  }

  @grok.decorators.panel({
    name: 'Biology | Admetica',
    meta: {role: 'widgets', domain: 'chem'},
  })
  static async admeticaWidget(
    @grok.decorators.param({ name: 'smiles', options: { semType: 'Molecule' }}) semValue: DG.SemanticValue): Promise<DG.Widget<any>> {
    const smiles = await grok.functions.call('Chem:convertMolNotation',
      { molecule: semValue.value, sourceNotation: DG.chem.Notation.Unknown, targetNotation: DG.chem.Notation.Smiles });
    return await getModelsSingle(smiles, semValue);
  }

  @grok.decorators.func()
  static async getModels(property?: string): Promise<string[]> {
    await setProperties();
    return properties.subgroup
      .filter((subg: Subgroup) => !property || subg.name === property)
      .flatMap((subg: Subgroup) => subg.models)
      .map((model: Model) => model.name);
  }

  @grok.decorators.func({ name: 'AdmeticaHT',
    meta: {role: 'hitTriageFunction'},
  })
  static async admeticaHT(
    table: DG.DataFrame,
    @grok.decorators.param({ options: { semType: 'Molecule' } }) molecules: DG.Column,
    @grok.decorators.param({ options: { choices: 'Admetica:getModels(\'Absorption\')', nullable: true } }) absorption: string[],
    @grok.decorators.param({ options: { choices: 'Admetica:getModels(\'Distribution\')', nullable: true } }) distribution: string[],
    @grok.decorators.param({ options: { choices: 'Admetica:getModels(\'Metabolism\')', nullable: true } }) metabolism: string[],
    @grok.decorators.param({ options: { choices: 'Admetica:getModels(\'Excretion\')', nullable: true } }) excretion: string[],
  ): Promise<void> {
    const models: string[] = [
      ...absorption,
      ...distribution,
      ...metabolism,
      ...excretion,
    ];
    await performChemicalPropertyPredictions(molecules, table, models);
  }

  @grok.decorators.editor({name: 'AdmeticaEditor'})
  static admeticaEditor(call: DG.FuncCall): void {
    const funcEditor = new AdmeticaBaseEditor();
    ui.dialog({ title: 'Admetica' })
      .add(funcEditor.getEditor())
      .onOK(async () => {
        const params = funcEditor.getParams();
        call.func
          .prepare({
            table: params.table,
            molecules: params.col,
            template: params.templateContent,
            models: params.models,
            addPiechart: params.addPiechart,
            addForm: params.addForm,
          })
          .call(true);
      })
      .show();
  }

  @grok.decorators.func({
    'name': 'AdmeticaMenu',
    'top-menu': 'Chem | Admetica | Сalculate...',
    'editor': 'Admetica:AdmeticaEditor',
  })
  static async admeticaMenu(
    @grok.decorators.param({options: { description: 'Input data table' }}) table: DG.DataFrame,
    @grok.decorators.param({options: { semType: 'Molecule' }}) molecules: DG.Column,
      template: string,
      models: string[],
      addPiechart: boolean,
      addForm: boolean,
  ): Promise<void> {
    await performChemicalPropertyPredictions(molecules, table, models, template, addPiechart, addForm);
  }

  @grok.decorators.func({
    name: 'getAdmeProperties',
    meta: {vectorFunc: 'true'},
    outputs: [{name: 'result', type: 'dataframe', options: {action: 'join(table)'}}],
  })
  static async getAdmeProperties(
    table: DG.DataFrame,
    @grok.decorators.param({ options: { semType: 'Molecule' } }) molecules: DG.Column,
    @grok.decorators.param({type: 'list<string>', options: { optional: true }}) props?: string[],
  ): Promise<DG.DataFrame> {
    const isMolblock = molecules.meta.units === DG.UNITS.Molecule.MOLBLOCK ||
      (!molecules.meta.units && DG.Detector.sampleCategories(molecules, (s) => s.includes('M  END'), 1));

    const values = new Array(molecules.length + 1);
    values[0] = molecules.name;
    for (let i = 0; i < molecules.length; i++) {
      const value = molecules.get(i);
      values[i + 1] = isMolblock ? `"${value}"` : value;
    }
    const csv = values.join('\n');

    // If no properties specified, use all available models
    const models = (props ?? await this.getModels()).join(',');
    let admeticaResults = await grok.functions.call('Admetica:run_admetica', {
      csv: csv,
      models: models,
      raiseException: false,
    });

    admeticaResults = await convertLD50(admeticaResults, molecules);
    return admeticaResults;
  }

  @grok.decorators.func({
    name: 'getAdmePropertiesSingle',
    description: 'Predicts ADME properties for a given molecule.',
  })
  static async getAdmePropertiesSingle(
    @grok.decorators.param({ options: { semType: 'Molecule' } }) molecule: string,
  ): Promise<DG.DataFrame> {
    const col = DG.Column.fromStrings('Molecule', [molecule]);
    return await PackageFunctions.getAdmeProperties(DG.DataFrame.fromColumns([col]), col);
  }

  @grok.decorators.app({name: 'Admetica', meta: {icon: 'images/vlaaivis.png', browsePath: 'Chem'}})
  static async runAdmeticaApplication(): Promise<DG.ViewBase | null> {
    const parent = grok.functions.getCurrentCall();
    return await initializeAdmeticaApp(true, parent);
  }

  @grok.decorators.demo({
    name: 'Admetica Demo',
    description: 'Evaluating ADMET properties',
    meta: {demoPath: 'Cheminformatics | Admetica'},
    outputs: [],
  })
  static async admeticaDemo(): Promise<DG.ViewBase | null> {
    return await initializeAdmeticaApp(false);
  }
}

async function initializeAdmeticaApp(addPath: boolean = true, parent: DG.FuncCall | null = null): Promise<DG.ViewBase | null> {
  const app = new AdmeticaViewApp(parent);
  app.addPath = addPath;
  await app.init();
  return grok.shell.addPreview(app.tableView!);
}
