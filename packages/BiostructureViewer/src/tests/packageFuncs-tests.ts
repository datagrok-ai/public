import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category/*, expect*/, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';

category('packageFuncs', () => {
  test('nglForGridTestApp', async () => {
    await grok.functions.call(`${_package.name}:nglForGridTestApp`, {});
  });

  test('nglViewerApp', async () => {
    await grok.functions.call(`${_package.name}:nglViewerApp`, {});
  });

  test('biostructureViewerApp', async () => {
    await grok.functions.call(`${_package.name}:biostructureViewerApp`, {});
  });

  test('biotrackViewerApp', async () => {
    await grok.functions.call(`${_package.name}:biotrackViewerApp`, {});
  });

  test('biostructureAndTrackViewerApp', async () => {
    await grok.functions.call(`${_package.name}:biostructureAndTrackViewerApp`, {});
  });
});
