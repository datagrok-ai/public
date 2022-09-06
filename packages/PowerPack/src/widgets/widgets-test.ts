import {after, before, category, expect, test} from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AboutWidget} from "./about-widget";
import {CommunityWidget} from "./community-widget";
import {HtmlWidget} from "./html-widget";
import {RecentProjectsWidget} from "./recent-projects-widget";
import {SystemStatusWidget} from "./system-status-widget";

category('Widgets', () => {

  test('AboutWidget', async () => new AboutWidget());
  test('CommunityWidget', async () => new CommunityWidget());
  test('HtmlWidget', async () => new HtmlWidget());
  test('RecentProjectsWidget', async () => new RecentProjectsWidget());
  test('SystemStatusWidget', async () => new SystemStatusWidget());

});
