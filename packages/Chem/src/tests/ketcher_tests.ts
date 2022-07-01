import {category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';
import {Func} from "datagrok-api/src/entities";
import * as grok from 'datagrok-api/grok';
const exampleSmiles = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3';
const convertedSmarts = '[#6]-[#6](-[#6](=O)-[#8]-[#6]-[#6]-[#6]-c1cccnc1)-c1cccc(c1)-[#6](=O)-c1ccccc1';

category('sketcher testing', () => {
    test('smiles-to-mol', async () => {
        const funcs = Func.find({tags: ['moleculeSketcher']});
        for (let f of funcs) {
            const smilesString = exampleSmiles;
            // @ts-ignore
            const sketcher = await f!.apply();
            sketcher.smiles = smilesString;
            const resultMol = sketcher.molFile;
            const convertedSmiles = grok.chem.convert(resultMol, 'mol', 'smiles');
            expect(convertedSmiles, smilesString);
        }

    });

    test('mol-to-smiles', async () => {
        const funcs = Func.find({tags: ['moleculeSketcher']});
        for (let f of funcs) {
            const smilesString = exampleSmiles;
            // @ts-ignore
            const sketcher = await f!.apply();
            const startMol = grok.chem.convert(smilesString, 'smiles', 'mol');
            sketcher.mol = startMol;
            const resultSmiles = sketcher.smiles;
            expect(resultSmiles, smilesString);
        }
    });
    test('mol-to-smarts', async () => {
        const funcs = Func.find({tags: ['moleculeSketcher']});
        for (let f of funcs) {
            const smilesString = exampleSmiles;
            // @ts-ignore
            const sketcher = await f!.apply();
            const startMol = grok.chem.convert(smilesString, 'smiles', 'mol');
            sketcher.mol = startMol;
            const resultSmarts = await sketcher.getSmarts();
            expect(resultSmarts, convertedSmarts);
        }
    });
    test('smiles-to-smarts', async () => {
        const funcs = Func.find({tags: ['moleculeSketcher']});
        for (let f of funcs) {
            const smilesString = exampleSmiles;
            // @ts-ignore
            const sketcher = await f!.apply();
            sketcher.smiles = smilesString;
            const resultSmarts = await sketcher.getSmarts();
            expect(resultSmarts, convertedSmarts);
        }

    });

});
