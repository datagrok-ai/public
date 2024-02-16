import {awaitCheck, before, category, delay, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { runAdmetox, performChemicalPropertyPredictions, addCalculationsToTable } from '../utils/admetox-utils';

category('Admetox', () => {
	let v: DG.TableView;
	let molecules: DG.DataFrame;
	let smilesColumn: DG.Column;
	
	before(async () => {
		grok.shell.closeAll();
		grok.shell.windows.showProperties = true;
  });
	
	test('Container', async () => {
		const admetDockerfile = await grok.dapi.docker.dockerContainers.filter('admetox').first();
		expect(admetDockerfile != null, true);
	});

	test('Container. Post request', async () => {
		const smiles = `smiles
		O=C1Nc2ccccc2C(C2CCCCC2)=NC1`;
		const bbbResults = await runAdmetox(smiles, 'PPB,VD,BBB', 'false');
		expect(bbbResults != null, true);
	})

  test('Calculate dialog. UI', async () => {
		molecules = grok.data.demo.molecules(100);
		v = grok.shell.addTableView(molecules);
		grok.shell.topMenu.find('Chem').group('ADME/Tox').find('Calculate...').click();
		await awaitCheck(() => document.querySelector('.d4-dialog') !== null,
			'cannot find ADME/Tox dialog', 5000);
		const admetoxDialog = document.querySelector('.d4-dialog');
		await awaitCheck(() => admetoxDialog!.querySelectorAll('.d4-tree-view-group-host > .d4-tree-view-group').length === 4,
			'properties number inside ADME/Tox dialog is different than expected', 5000);
		const smilesColumn = admetoxDialog!.querySelector('.d4-column-selector-column') as HTMLElement;
		await awaitCheck(() => smilesColumn!.innerText === 'smiles',
		  'column inside ADME/Tox dialog is different than expected', 5000);
		await delay(1000);
		const expandBtn = admetoxDialog!.querySelectorAll('.d4-tree-view-tri')[0] as HTMLElement;
		expandBtn.click();
		await delay(2000);
		expectArray(Array.from(admetoxDialog!.querySelectorAll('.d4-tree-view-item-label')).map((item) => item.innerHTML), 
		  ['Pgp-Inhibitor', 'Pgp-Substrate', 'HIA', 'F(20%)', 'F(30%)', 'BBB', 'CYP1A2-Inhibitor', 'CYP1A2-Substrate',
		  'CYP3A4-Inhibitor', 'CYP3A4-Substrate', 'CYP2C19-Inhibitor', 'CYP2C19-Substrate', 'CYP2C9-Inhibitor', 'CYP2C9-Substrate',
		  'CYP2D6-Inhibitor', 'CYP2D6-Substrate', 'Ames', 'SkinSen']);
		v.close();
		grok.shell.o = ui.div();
	});

	test('Calculate dialog. Added properties', async () => {
		molecules = grok.data.demo.molecules(10);
		v = grok.shell.addTableView(molecules);
		smilesColumn = molecules.columns.bySemType(DG.SEMTYPE.MOLECULE)!;
		const newTableColumn = 'Pgp-Inhibitor';
		const predictions = addCalculationsToTable(molecules);
		await delay(1000);
		const admetoxDialog = document.querySelector('.d4-dialog');
		(admetoxDialog!.querySelectorAll('input[type="checkbox"]')[0] as HTMLElement).click();
		await delay(1000);
		(admetoxDialog!.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement).click();	
		await predictions;
		expect(molecules.columns.names().includes(newTableColumn), true, `${newTableColumn} column has not been added`);
		expect(molecules.col(newTableColumn)!.get(0), 0.44976675510406494, `Calculated value for ${newTableColumn} is incorrect`);
		expect(molecules.col(newTableColumn)!.colors.getColor(0), 4293426297, 'Wrong color coding was added');
		expect(molecules.col(newTableColumn)!.colors.getColor(4), 4282627449, 'Wrong color coding was added');
	});

	test('Calculate. For single cell', async () => {
		molecules = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(molecules);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;
    await awaitPanel(pp, 'ADME/Tox', 3000);
    (document.querySelector('.fa-chevron-square-down') as HTMLElement)?.click();
    const absorp = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Absorption') as HTMLElement;
    if (!absorp.classList.contains('expanded')) {
			absorp.click();
		}
		const absorpRes = `Pgp-Inhibitor\t0.45
		Pgp-Substrate\t0.50
		HIA\t0.64
		F(20%)\t0.57
		F(30%)\t0.88`;
		await awaitCheck(() => (pp.getElementsByClassName('d4-table d4-item-table d4-info-table')[3] as HTMLElement).innerText === absorpRes, 'Results for single cell differ', 8000);
	});

	test('Calculate.Benchmark', async () => {
		molecules = DG.Test.isInBenchmark 
		  ? await grok.data.files.openTable("Demo:Files/chem/smiles_1M.zip")
			: grok.data.demo.molecules();
		DG.time('ADME/Tox post request', () => runAdmetox(molecules.toCsv(), 'Pgp-Inhibitor,Pgp-Substrate,HIA,F(20%),F(30%)', 'false'));
	})
});

async function awaitPanel(pp: HTMLElement, name: string, ms: number = 5000): Promise<void> {
  await awaitCheck(() => {
    return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === name) !== undefined;
  }, `cannot find ${name} property`, ms);
}