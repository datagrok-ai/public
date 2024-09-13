import {awaitCheck, before, category, delay, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { runAdmetica, performChemicalPropertyPredictions, getQueryParams, properties, setProperties } from '../utils/admetica-utils';
import { fetchWrapper } from '@datagrok-libraries/utils/src/fetch-utils';

category('Admetica', () => {
  let v: DG.TableView;
  let molecules: DG.DataFrame;
  let smilesColumn: DG.Column;
    
  before(async () => {
    grok.shell.closeAll();
    grok.shell.windows.showProperties = false;
    await setProperties();
  });
    
  test('Container', async () => {
    const admetDockerfile = await grok.dapi.docker.dockerContainers.filter('admetica').first();
    expect(admetDockerfile != null, true);
  });

  test('Container. Post request', async () => {
    const smiles = `smiles
    O=C1Nc2ccccc2C(C2CCCCC2)=NC1`;
    const bbbResults = await fetchWrapper(() => runAdmetica(smiles, 'PPBR,VDss', 'false'));
    expect(bbbResults != null, true);
  }, {timeout: 100000});

  test('Calculate dialog. UI', async () => {
    molecules = grok.data.demo.molecules(100);
    v = grok.shell.addTableView(molecules);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    grok.shell.topMenu.find('Chem').group('Admetica').find('Сalculate...').click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length > 0, 'cannot open Admetica dialog', 2000);
    const admeticaDialog = returnDialog('Admetica')?.root;
    const settingsIcon = admeticaDialog?.querySelector('.grok-icon.grok-font-icon-settings') as HTMLElement;
    settingsIcon.click();
    await awaitCheck(() => admeticaDialog!.querySelectorAll('.d4-tree-view-group-host > .d4-tree-view-group').length === 4,
      'properties number inside Admetica dialog is different than expected', 5000);
    const smilesColumn = admeticaDialog!.querySelector('.d4-column-selector-column') as HTMLElement;
    await awaitCheck(() => smilesColumn!.innerText === 'smiles',
      'column inside Admetica dialog is different than expected', 5000);
    const models = properties.subgroup
      .flatMap((subgroup: any) => subgroup.models
      .map((model: any) => model.name));
    expectArray(Array.from(admeticaDialog!.querySelectorAll('.d4-tree-view-item-label')).map((item) => item.innerHTML), models);
    v.close();
    grok.shell.o = ui.div();
  });

  test('Calculate dialog. Added properties', async () => {
    molecules = grok.data.demo.molecules(5);
    v = grok.shell.addTableView(molecules);
    await delay(1000);
    smilesColumn = molecules.columns.bySemType(DG.SEMTYPE.MOLECULE)!;
    const newTableColumn = 'Caco2';
    await performChemicalPropertyPredictions(smilesColumn, v.dataFrame, newTableColumn);
    expect(molecules.columns.names().includes(newTableColumn), true, `${newTableColumn} column has not been added`);
    expect(molecules.col(newTableColumn)!.get(0), -4.615971565246582, `Calculated value for ${newTableColumn} is incorrect`);
    expect(molecules.col(newTableColumn)!.meta.colors.getColor(0), 4278255360, 'Wrong color coding was added');
    expect(molecules.col(newTableColumn)!.meta.colors.getColor(4), 4288177664, 'Wrong color coding was added');
  }, {timeout: 100000});

  test('Calculate. For single cell', async () => {
    molecules = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(molecules);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    grok.shell.windows.showProperties = true;
    const table = v.dataFrame;
    table.currentCell = table.cell(0, 'smiles');
    await delay(1000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;
    await awaitPanel(pp, 'Admetica', 6000);
    (document.querySelector('.fa-chevron-square-down') as HTMLElement)?.click();
    const distribution = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Distribution') as HTMLElement;
    if (!distribution.classList.contains('expanded')) {
      distribution.click();
    }
    await delay(1000);
    const distributionRes = `
      PPBR\t82.26
      VDss\t8.48`;
    await awaitCheck(() => (pp.getElementsByClassName('d4-table d4-item-table d4-info-table')[2] as HTMLElement).innerText === distributionRes, 'Results for single cell differ', 8000);
  }, {timeout: 100000});

  test('Calculate.Benchmark column', async () => {
    const runAdmeticaBenchmark = async (moleculesCount: number) => {
        const molecules = grok.data.demo.molecules(moleculesCount);
        molecules.columns.remove('logD');
        const iterations = DG.Test.isInBenchmark ? 100 : 5;
        const args = [molecules.toCsv(), await getQueryParams(), 'false'];
        return await runInLoop(iterations, runAdmetica, ...args);
    };

    const mol1k = await runAdmeticaBenchmark(1000);
    const mol5k = await runAdmeticaBenchmark(5000);
    const mol10k = await runAdmeticaBenchmark(10000);

    return DG.toDart({"1k molecules": mol1k, "5k molecules": mol5k, "10k molecules": mol10k});
}, {timeout: 10000000000, benchmark: true });

  test('Calculate.Benchmark cell', async () => {
    const smiles = `smiles
    O=C1Nc2ccccc2C(C2CCCCC2)=NC1`;
    const iterations = DG.Test.isInBenchmark ? 100 : 10;
    const distributionSubgroup = properties.subgroup.find((subgroup: any) => subgroup.name === "Distribution");
    const distributionModels = distributionSubgroup ? distributionSubgroup.models.map((model: any) => model.name) : [];
    const args = [smiles, distributionModels, 'false'];
    const cellResults = await runInLoop(iterations, runAdmetica, ...args);
    return DG.toDart({"results": cellResults});
  }, {timeout: 1000000, benchmark: true});
});

async function runInLoop(iterations: number, func: (...args: string[]) => Promise<string | null>, ...args: string[]) {
  const results = new Array<number>(iterations);
  for (let i = 0; i < iterations; ++i) {
    const startTime = performance.now();
    await func(...args);
    const endTime = performance.now();
    results[i] = (endTime - startTime) / 1000;
  }
  const sum = results.reduce((p, c) => p + c, 0);
  return {'Iterations' : results.length, 'Average time': sum / results.length,
    'Min time': Math.min(...results), 'Max time': Math.max(...results)};
}

async function awaitPanel(pp: HTMLElement, name: string, ms: number = 5000): Promise<void> {
  await awaitCheck(() => {
    return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === name) !== undefined;
  }, `cannot find ${name} property`, ms);
}

export function returnDialog(dialogTitle: string): DG.Dialog | undefined {
  let dialog: DG.Dialog | undefined;
  for (let i = 0; i < DG.Dialog.getOpenDialogs().length; i++) {
    if (DG.Dialog.getOpenDialogs()[i].title == dialogTitle) {
      dialog = DG.Dialog.getOpenDialogs()[i];
      return dialog;
    }
  }
}