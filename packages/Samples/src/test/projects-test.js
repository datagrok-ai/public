var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
import * as grok from 'datagrok-api/grok';
import { category, test } from '@datagrok-libraries/utils/src/test';
category('Projects', () => {
    test('coffee_company', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('coffee_company');
    }));
    test('chem_demo', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('chem_demo');
    }), { skipReason: 'GROK-13612' });
    test('chem_tsne', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('chem_tsne');
    }), { skipReason: 'GROK-13621' });
    test('demo_notebooks', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('demo_notebooks');
    }), { skipReason: 'GROK-15481' });
    test('demog', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('demog');
    }));
    // test('ecg', async () => {
    //     await grok.dapi.projects.open('ecg')
    // });
    test('eeg', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('eeg');
    }));
    test('game-of-thrones', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('game_of_thrones');
    }), { skipReason: 'GROK-16289' });
    test('models', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('models');
    }));
    // test('pa-income-by-county', async () => {
    //     await grok.dapi.projects.open('pa_income_by_county')
    // });
    test('it-network', () => __awaiter(void 0, void 0, void 0, function* () {
        yield grok.dapi.projects.open('it_network');
    }));
    // test('plates', async () => {
    //     await grok.dapi.projects.open('plates')
    // });
    // test('time-series-decomposition', async () => {
    //     await grok.dapi.projects.open('time_series_decomposition')
    // });
});
//# sourceMappingURL=projects-test.js.map