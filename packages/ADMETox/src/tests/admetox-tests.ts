import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { accessServer } from '../admet-analysis/admet-calculation';
import { models } from '../admet-analysis/const';

category('Admetox', () => {
    let dockerId: string;
    
    before(async () => {
        dockerId = (await grok.dapi.dockerfiles.filter('admetox').first()).id;
    });
    
    test('dockerfile work', async () => {
        expect(await grok.dapi.dockerfiles.stop(dockerId), true);
        expect(await grok.dapi.dockerfiles.run(dockerId), true);
    });

    test('predict results', async () => {
        const testSmile = `smiles
        O=C1Nc2ccccc2C(C2CCCCC2)=NC1`;
        const csvResult = await accessServer(testSmile, Object.keys(models).toString());
        console.log(csvResult);
    })
});