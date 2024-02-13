import * as DG from 'datagrok-api/dg';
import { DockingPackage } from '../package-utils';
import { AutoDockDataType } from '../apps/auto-dock-app';

export const _package = new DockingPackage();
export const TARGET_PATH = 'System:AppData/Docking/targets';
export const CACHED_DOCKING: DG.LruCache<AutoDockDataType, DG.DataFrame> = new DG.LruCache<AutoDockDataType, DG.DataFrame>();
export const CACHED_MOLSTAR: DG.LruCache<string, DG.Widget> = new DG.LruCache<string, DG.Widget>();
export const AFFINITY_COL = 'affinity';
export const POSE_COL = 'pose';