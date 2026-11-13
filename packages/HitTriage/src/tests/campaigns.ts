import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package';
import {loadCampaigns} from '../app/utils';
import {AppName} from '../app/types';
import {CampaignJsonName} from '../app/consts';
import {defaultPermissions} from '../app/shared/permissions-editor';

category('Campaigns', () => {
  test('List skips unreadable campaigns and shows those shared with All users', async () => {
    const appName = `Tests ${DG.Utils.randomString(8)}` as AppName;
    const otherAuthorId = '00000000-0000-0000-0000-000000000001';
    try {
      await _package.files.writeAsText(`${appName}/campaigns/GOOD-1/${CampaignJsonName}`,
        JSON.stringify({name: 'GOOD-1', authorUserId: otherAuthorId, permissions: defaultPermissions}));
      await _package.files.writeAsText(`${appName}/campaigns/BAD-1/${CampaignJsonName}`, '{not json');
      const campaigns = await loadCampaigns(appName, []);
      expect(Object.keys(campaigns).join(','), 'GOOD-1');
    } finally {
      await _package.files.delete(appName);
    }
  });
});
