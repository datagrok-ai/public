/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { welcomeView } from "./welcome-view";
import { compareColumns } from './compare-columns';
import { DistributionProfilerViewer } from './distribution-profiler';
import {SystemStatusWidget} from "./widgets/system-status-widget";
import {RecentProjectsWidget} from "./widgets/recent-projects-widget";

export let _package = new DG.Package();

//name: test
//tags: autostart
export function test(): void {
  welcomeView();
}

//name: compareColumns
//top-menu: Data | Compare Columns...
export function _compareColumns(): void {
  compareColumns();
}

//name: distributionProfiler
//tags: viewer
//output: viewer result
export function _distributionProfiler(): DistributionProfilerViewer {
  return new DistributionProfilerViewer();
}

//name: welcomeView
//tags: autostart
export function _welcomeView(): void {
  welcomeView();
}

//output: widget result
export function systemStatusWidget(): DG.Widget {
  return new SystemStatusWidget();
}

//output: widget result
export function recentProjectsWidget(): DG.Widget {
  return new RecentProjectsWidget();
}
