import { category, expect, expectObject, test } from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';

category('UserSettingsStorage', () => {
  test('Add and get value', async () => {
    const storageName = 'js-api-storage-name1';
    const key = 'postValueKey';
    const value = 'value';
    grok.userSettings.add(storageName, key, value);
    expect(grok.userSettings.getValue(storageName, key), value);
  });

  test('AddAll and get map', async () => {
    const storageName = 'js-api-storage-name2';
    const value = { 'key': 'value' };
    grok.userSettings.addAll(storageName, value);
    const receivedValue = grok.userSettings.get(storageName);
    expect(JSON.stringify(receivedValue), JSON.stringify(value));
  });

  test('put', async () => {
    const storageName = 'js-api-storage-name3';
    const value1 = { 'key1': 'value1' };
    const value2 = { 'key2': 'value2' };

    grok.userSettings.addAll(storageName, value1);
    grok.userSettings.addAll(storageName, value2);
    expectObject(grok.userSettings.get(storageName) ?? {}, { 'key1': 'value1', 'key2': 'value2' });
    grok.userSettings.put(storageName, value2);
    expectObject(grok.userSettings.get(storageName) ?? {}, value2);
  });

  test('delete', async () => {
    const storageName = 'js-api-storage-name4';
    const key = 'postValueKey';
    const value = 'value';

    grok.userSettings.add(storageName, key, value);
    grok.userSettings.delete(storageName, key);

    const receivedValue = grok.userSettings.getValue(storageName, key);
    expect(receivedValue == undefined);
  }, { stressTest: true });

  test('credentials', async () => {
    console.log(_package);
    let url = `https://${window.location.host}/api/credentials/for/${_package.name}`;
    const headers = new Headers();
    headers.append('Accept', 'application/json');
    headers.append('testAuth', `test`);
    headers.append('original-url', url);
    headers.append('original-method', 'POST');

    const requestOptions: RequestInit = {
      method: "POST",
      headers: headers,
      redirect: "follow",

      body: JSON.stringify({test: 'test'}),
    };
    const response = await grok.dapi.fetchProxy(url, requestOptions);
    
    let credentialsParams = (await _package.getCredentials()).parameters;

    expect(credentialsParams['test'], 'test');
  }, {skipReason: 'Skipped for 1.25.0'});

}, { owner: 'ppolovyi@datagrok.ai' });
