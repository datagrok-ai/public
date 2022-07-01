import {category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';
import {Func} from "datagrok-api/src/entities";
import * as grok from 'datagrok-api/grok';


category('sketcher testing', () => {
    test('smiles-to-mol', async () => {
        const funcs = Func.find({tags: ['moleculeSketcher']});
        for (let f of funcs) {
            const smilesString = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3'
            // @ts-ignore
            const sketcher = await f!.apply();
            sketcher.smiles = smilesString;
            const resultMol = sketcher.molFile;
            const convertedSmiles = grok.chem.convert(resultMol, 'mol', 'smiles');
            expect(convertedSmiles, smilesString);
        }

    });

    test('mol-to-smiles', async () => {
     
    });

});
