import {after, before, category, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi: user data storage', () => {

    test('Dapi: user data storage - post and get value', async () => {
        let storageName = 'js-api-storage-name1';
        let key = 'postValueKey';
        let value = 'value';
        await grok.dapi.userDataStorage.postValue(storageName, key, value);
        expect(await grok.dapi.userDataStorage.getValue(storageName, key), value);
    });

    test('Dapi: user data storage - post and get map', async () => {
        let storageName = 'js-api-storage-name2';
        let value = {'key': 'value'};
        await grok.dapi.userDataStorage.post(storageName, value);
        let receivedValue = await grok.dapi.userDataStorage.get(storageName);
        expect(JSON.stringify(receivedValue), JSON.stringify(value));
    });

    test('Dapi: user data storage - put', async () => {
        let storageName = 'js-api-storage-name3';
        let value1 = {'key': 'value1'};
        let value2 = {'key': 'value2'};

        await grok.dapi.userDataStorage.post(storageName, value1);
        await grok.dapi.userDataStorage.put(storageName, value2);

        let receivedValue = await grok.dapi.userDataStorage.get(storageName);
        expect(JSON.stringify(receivedValue), JSON.stringify(value2));
    });

    test('Dapi: user data storage - put', async () => {
        let storageName = 'js-api-storage-name4';
        let key = 'postValueKey';
        let value = 'value';

        await grok.dapi.userDataStorage.postValue(storageName, key, value);
        await grok.dapi.userDataStorage.remove(storageName, key);

        let receivedValue = await grok.dapi.userDataStorage.getValue(storageName, key);
        expect(receivedValue, '');
    });

});
