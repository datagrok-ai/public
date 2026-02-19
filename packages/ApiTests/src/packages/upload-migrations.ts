import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {before, category, expect, test,  expectExceptionAsync, after} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';

category('Packages: migrations', () => {

  // we want to have the migration passed
  let packageWithMigration: Uint8Array;
  // we want to have the migration passed
  let packageWithSecondMigration: Uint8Array;
  // we want to throw an error on permission
  let packageWithDropEntities: Uint8Array;
  let key: string;

  before(async () => {
    packageWithMigration = await _package.files.readAsBytes('package/packageWithMigration.zip');
    packageWithSecondMigration = await _package.files.readAsBytes('package/packageWithSecondMigration.zip');
    packageWithDropEntities = await _package.files.readAsBytes('package/packageWithDropEntities.zip');
    const response = await fetch(`${grok.dapi.root}/users/current/dev_key?getNew=false`,
      {method: 'GET', credentials: 'include'});
    key = await response.json();

    await publish(packageWithMigration, 'apitestsdb', false);    
  });

  after(async () => {
    await deletePackage('apitestsdb');
  });

  async function publish(packageData: Uint8Array, name: string, debug: boolean) {
    const uploadResponse = await fetch(`${grok.dapi.root}/packages/dev/${key}/${name}?debug=${debug}&rebuild=false`,
      {method: 'POST', body: packageData as BodyInit});
    expect(uploadResponse.status, 200);
    let text = await uploadResponse.text();
    expect(text.indexOf('ApiError'), -1);
    console.log('wait 10 seconds to fully init package');
    await DG.delay(10000);
  }

  async function deletePackage(packageName: string) {
    const uploadResponse = await fetch(`${grok.dapi.root}/packages/dev/${key}/${packageName}`,
      {method: 'DELETE'});
    let text = await uploadResponse.text();
    expect(text.indexOf('ApiError'), -1);
  }

  test('Query uploaded package', async () => {
    await grok.functions.call('apitestsdb:AccessTableA');
  });

  test('Upload new migration', async () => {
    // in release: published, updated
    await publish(packageWithSecondMigration, 'apitestsdb', false);
    await grok.functions.call('apitestsdb:AccessTableB');
    
    // previous versions should be disabled 
    let versionsList = await grok.dapi.packages.filter('shortName = "Apitestsdb"').list();
    expect(versionsList.length, 1);
  }, { timeout: 50000 });

  test('Isolation', async () => {
    await expectExceptionAsync(() => publish(packageWithDropEntities, 'apitestsbad', false));
    await deletePackage('apitestsbad');
  });
}, {
  owner: 'aparamonov@datagrok.ai'
});
