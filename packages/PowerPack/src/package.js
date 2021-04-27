import * as DG from 'datagrok-api/dg';
import { welcomeView } from "./welcome-view";
import { compareColumns } from './compare-columns';
export let _package = new DG.Package();
//name: test
//tags: autostart
export function test() {
    welcomeView();
}
//name: compareColumns
//top-menu: Data | Compare Columns
export function _compareColumns() {
    compareColumns();
}
