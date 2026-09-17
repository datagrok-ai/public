import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {AboutWidget} from '../widgets/about-widget';
import {CommunityWidget} from '../widgets/community-widget';
import {HtmlWidget} from '../widgets/html-widget';
import {renderPlaylist} from '../widgets/learning-widget';
import {RecentProjectsWidget} from '../widgets/recent-projects-widget';

category('Widgets', () => {
  test('AboutWidget', async () => {new AboutWidget();});
  test('CommunityWidget', async () => {new CommunityWidget();});
  test('HtmlWidget', async () => {new HtmlWidget();});
  test('RecentProjectsWidget', async () => {new RecentProjectsWidget();});
  test('LearningWidget play icon opens playlist', async () => {
    const row = renderPlaylist({id: 'PL123', title: 'Meetings', description: 'd', color: '#c9b6b6'});
    const open = window.open;
    let opened = '';
    window.open = (url?: string | URL) => {opened = String(url); return null;};
    try {
      (row.querySelector('.power-pack-learn-icon') as HTMLElement).click();
    } finally {
      window.open = open;
    }
    expect(opened, 'https://www.youtube.com/playlist?list=PL123');
  });
});
