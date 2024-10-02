import * as grok from 'datagrok-api/grok';
import {category, test} from '@datagrok-libraries/utils/src/test';


category('Projects', () => {
    test('coffee_company', async () => {
        await grok.dapi.projects.open('coffee_company')
    });
    test('chem_demo', async () => {
        await grok.dapi.projects.open('chem_demo')
    }, {skipReason: 'GROK-13612'});
    test('chem_tsne', async () => {
        await grok.dapi.projects.open('chem_tsne')
    }, {skipReason: 'GROK-13621'});
    test('demo_notebooks', async () => {
        await grok.dapi.projects.open('demo_notebooks')
    }, {skipReason: 'GROK-15481'});
    test('demog', async () => {
        await grok.dapi.projects.open('demog')
    });
    // test('ecg', async () => {
    //     await grok.dapi.projects.open('ecg')
    // });
    test('eeg', async () => {
        await grok.dapi.projects.open('eeg')
    });
    test('game-of-thrones', async () => {
        await grok.dapi.projects.open('game_of_thrones')
    }, {skipReason: 'GROK-16289'});
    test('models', async () => {
        await grok.dapi.projects.open('models')
    });
    // test('pa-income-by-county', async () => {
    //     await grok.dapi.projects.open('pa_income_by_county')
    // });
    test('it-network', async () => {
        await grok.dapi.projects.open('it_network')
    });
    // test('plates', async () => {
    //     await grok.dapi.projects.open('plates')
    // });
    // test('time-series-decomposition', async () => {
    //     await grok.dapi.projects.open('time_series_decomposition')
    // });
});
