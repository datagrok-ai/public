import {category, expect, expectFloat, test, delay} from '@datagrok-libraries/utils/src/test';
import {Func} from "datagrok-api/src/entities";
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {chem} from 'datagrok-api/grok';
import Sketcher = chem.Sketcher;
import { _package } from '../package-test';
const exampleSmiles = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3';
const convertedSmarts = '[#6]-[#6](-[#6](=O)-[#8]-[#6]-[#6]-[#6]-c1cccnc1)-c1cccc(c1)-[#6](=O)-c1ccccc1';

category('sketcher testing', () => {
    test('smiles-to-mol', async () => {
        const module = await grok.functions.call('Chem:getRdKitModule');
        const funcs = Func.find({tags: ['moleculeSketcher']});
        const mol = module.get_mol(exampleSmiles);
        const sw = new Sketcher();
        const dg = ui.dialog().add(sw).show();
        await initSketcher(sw);
        for (let func of funcs) {         
            const fn = func.friendlyName;
            await sw.setSketcher(fn, exampleSmiles);
            const t = new Promise((resolve, reject) => {
                sw.onChanged.subscribe(async (_: any) => {
                    try {
                        const resultMol = sw.getMolFile();
                        resolve(resultMol);
                    } catch (error) {
                        reject(error);
                    }
                });
              });
            const resMolblock = await t;
            const mol2 = module.get_mol(resMolblock);
            const match1 = mol.get_substruct_match(mol2);
            expect(match1 !== '{}', true);
            const match2 = mol2.get_substruct_match(mol);
            expect(match2 !== '{}', true);
            mol2.delete();
        }
        mol?.delete();
        dg.close();
    });

    /*test('mol-to-smarts', async () => {
        const funcs = Func.find({tags: ['moleculeSketcher']});
        const mol = (await grok.functions.call('Chem:getRdKitModule')).get_qmol(convertedSmarts);

        //getMolblock
        const sw = new Sketcher();
        const dg = ui.dialog().add(sw).show();
        await initSketcher(sw);
        for (let func of funcs) {         
            const fn = func.friendlyName;
            await sw.setSketcher(fn, exampleSmiles);
            const molblock = mol.get_molblock();
            sw.setMolFile(molblock);
            const t = new Promise((resolve, reject) => {
                sw.onChanged.subscribe(async (_: any) => {
                    try {
                        const resultSmarts = sw.getSmarts();
                        resolve(resultSmarts);
                    } catch (error) {
                        reject(error);
                    }
                });
              });
            const resSmarts = await t;
            expect(resSmarts, convertedSmarts);
        }
        mol?.delete();
        dg.close();
         
    });*/

    /*test('mol-to-smiles', async () => { 
        const module = await grok.functions.call('Chem:getRdKitModule');
        const funcs = Func.find({tags: ['moleculeSketcher']});
        const mol = module.get_mol(exampleSmiles);
        const sw = new Sketcher();
        const dg = ui.dialog().add(sw).show();
        await initSketcher(sw);
        for (let func of funcs) {         
            const fn = func.friendlyName;
            const molblock = mol.get_molblock()
            await sw.setSketcher(fn, molblock);
            const t = new Promise((resolve, reject) => {
                sw.onChanged.subscribe(async (_: any) => {
                    try {
                        const resultSmiles = sw.getSmiles();
                        resolve(resultSmiles);
                    } catch (error) {
                        reject(error);
                    }
                });
              });
            const resSmiles = await t;
            const mol2 = module.get_mol(resSmiles);
            expect(mol.get_smiles(), mol2.get_smiles());
            mol2.delete();
        }
        mol?.delete();
        dg.close();
    });

    test('smiles-to-smarts', async () => {
        const funcs = Func.find({tags: ['moleculeSketcher']});
        const sw = new Sketcher();
        const dg = ui.dialog().add(sw).show();
        await initSketcher(sw);
        for (let func of funcs) {         
            const fn = func.friendlyName;
            await sw.setSketcher(fn, exampleSmiles);
            const t = new Promise((resolve, reject) => {
                sw.onChanged.subscribe(async (_: any) => {
                    try {
                        const resultSmarts = sw.getSmarts();
                        resolve(resultSmarts);
                    } catch (error) {
                        reject(error);
                    }
                });
              });
            const resSmart = await t;
            expect(resSmart, convertedSmarts);
        }
        dg.close();
    });
    
    test('molfileV2000', async () => {
        const data = DG.DataFrame.fromCsv(await _package.files.readAsText('test.csv'));
        const molfileV2000 = data.get('molecule', 0);
        const funcs = Func.find({tags: ['moleculeSketcher']});
        const sw = new Sketcher();
        const dg = ui.dialog().add(sw).show();
        await initSketcher(sw);
        for (let func of funcs) {         
            const fn = func.friendlyName;
            await sw.setMolFile(molfileV2000);
            const t = new Promise((resolve, reject) => {
                sw.onChanged.subscribe(async (_: any) => {
                    try {
                        const resultMolfile = sw.getMolFile();
                        resolve(resultMolfile);
                    } catch (error) {
                        reject(error);
                    }
                });
              });
            let resMolfile = await t;
            expect(resMolfile, molfileV2000);
        }
        dg.close();
    });
    
    test('smarts', async () => {
        const data = DG.DataFrame.fromCsv(await _package.files.readAsText('test-consistency-smarts-mol.csv'));
        const smarts = data.get('SMARTS', 0);
        const funcs = Func.find({tags: ['moleculeSketcher']});
        const sw = new Sketcher();
        const dg = ui.dialog().add(sw).show();
        await initSketcher(sw);
        for (let func of funcs) {         
            //const fn = func.friendlyName;
            await sw.setSmarts(smarts);
            const t = new Promise((resolve, reject) => {
                sw.onChanged.subscribe(async (_: any) => {
                    try {
                        const resultSmarts = sw.getSmarts();
                        resolve(resultSmarts);
                    } catch (error) {
                        reject(error);
                    }
                });
              });
            const resSmart = await t;
            expect(resSmart, smarts);
        }
        dg.close();
    });
    
    test('molfileV3000', async () => {
        const data = DG.DataFrame.fromCsv(await _package.files.readAsText('v3000_sample.csv'));
        const molfileV3000 = data.get('molecule', 0);
        const funcs = Func.find({tags: ['moleculeSketcher']});
        const sw = new Sketcher();
        const dg = ui.dialog().add(sw).show();
        await initSketcher(sw);
        for (let func of funcs) {         
            const fn = func.friendlyName;
            await sw.setMolFile(molfileV3000);
            const t = new Promise((resolve, reject) => {
                sw.onChanged.subscribe(async (_: any) => {
                    try {
                        const resultMolfile = sw.getMolFile();
                        resolve(resultMolfile);
                    } catch (error) {
                        reject(error);
                    }
                });
              });
            let resMolfile = await t;
            expect(resMolfile, molfileV3000);
        }
        dg.close();
    });*/

});

async function initSketcher(sw: Sketcher) {
    const t = new Promise((resolve, reject) => {
        sw.onChanged.subscribe(async (_: any) => {
            try {
                resolve(true);
            } catch (error) {
                reject(error);
            }
        });
      });
    sw.setSmiles(exampleSmiles);
    await t;
}
