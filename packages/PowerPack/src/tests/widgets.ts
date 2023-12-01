import {category, test} from '@datagrok-libraries/utils/src/test';
import {AboutWidget} from '../widgets/about-widget';
import {CommunityWidget} from '../widgets/community-widget';
import {HtmlWidget} from '../widgets/html-widget';
import {RecentProjectsWidget} from '../widgets/recent-projects-widget';
import {SystemStatusWidget} from '../widgets/system-status-widget';

category('Widgets', () => {
  test('AboutWidget', async () => {new AboutWidget();});
  test('CommunityWidget', async () => {new CommunityWidget();});
  test('HtmlWidget', async () => {new HtmlWidget();});
  test('RecentProjectsWidget', async () => {new RecentProjectsWidget();});
  test('SystemStatusWidget', async () => {new SystemStatusWidget();});
}, {clear: false});
