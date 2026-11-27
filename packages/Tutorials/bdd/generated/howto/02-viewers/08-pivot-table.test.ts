/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/08-pivot-table.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {pickInOpenColumnPicker} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clearSavedParameters, pivotedAggregationMatches} from '@datagrok-libraries/bdd/bindings/tiers/viewers/pivot-table';
import {addViewer, clickArea, closeContextMenu, noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragAreaOntoWidget} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/08-pivot-table.feature", import.meta.url);
  test("Summarize a table with a pivot table", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await session.step(23, "And user adds a pivot table viewer", () => addViewer(page, "pivot table"));
    await session.step(24, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(25, "And the \"pivot\" reading of pivot table viewer should be \"SEVERITY\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEVERITY"));
    await session.step(26, "When the video says \"A pivot table sums up the rows of a table by group\"", () => videoSays(page, "A pivot table sums up the rows of a table by group"));
    await session.step(28, "And user clicks on the \"add group by\" area of pivot table viewer", () => clickArea(page, "add group by", el("pivot table viewer")), undefined, "Group by: one row per value");
    await session.step(29, "And user picks \"RACE\" in the open column picker", () => pickInOpenColumnPicker(page, "RACE"));
    await session.step(30, "Then the \"group by\" reading of pivot table viewer should be \"DIS_POP, RACE\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP, RACE"));
    await session.step(32, "When user clicks on the \"remove group by chip DIS_POP\" area of pivot table viewer", () => clickArea(page, "remove group by chip DIS_POP", el("pivot table viewer")), undefined, "The x on a chip takes it out");
    await session.step(33, "Then the \"group by\" reading of pivot table viewer should be \"RACE\"", () => readingReads(page, "group by", el("pivot table viewer"), "RACE"));
    await session.step(34, "And the \"aggregated rows\" reading of pivot table viewer should be 4", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 4));
    await session.step(36, "When user clicks on the \"remove pivot chip SEVERITY\" area of pivot table viewer", () => clickArea(page, "remove pivot chip SEVERITY", el("pivot table viewer")), undefined, "Pivot: one column per value");
    await session.step(37, "And user clicks on the \"add pivot\" area of pivot table viewer", () => clickArea(page, "add pivot", el("pivot table viewer")));
    await session.step(38, "And user picks \"SEX\" in the open column picker", () => pickInOpenColumnPicker(page, "SEX"));
    await session.step(39, "Then the \"pivot\" reading of pivot table viewer should be \"SEX\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEX"));
    await session.step(40, "And the aggregated values of pivot table viewer should match \"avg(AGE)\" grouped by \"RACE\" pivoted on \"SEX\"", () => pivotedAggregationMatches(page, el("pivot table viewer"), "avg(AGE)", "RACE", "SEX"));
    await session.step(42, "When user drags the \"header HEIGHT\" area of grid onto the \"aggregate row\" area of pivot table viewer", () => dragAreaOntoWidget(page, "header HEIGHT", el("grid"), "aggregate row", el("pivot table viewer")), undefined, "Aggregate: what to compute. Or drag a column from the grid");
    await session.step(43, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(HEIGHT)"));
    await session.step(45, "When user picks \"Aggregation > max\" from the context menu of the \"aggregate chip avg(AGE)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > max", "aggregate chip avg(AGE)", el("pivot table viewer")), undefined, "Right-click a chip to change its function");
    await session.step(46, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(47, "Then the \"aggregate\" reading of pivot table viewer should be \"max(AGE), avg(HEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "max(AGE), avg(HEIGHT)"));
    await session.step(49, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")), undefined, "ADD puts the result into a new table");
    await session.step(50, "Then table \"demog aggregation\" should be open", () => tableOpen(page, "demog aggregation"));
    await session.step(51, "And the \"demog aggregation\" view should be current", () => viewIsCurrent(page, "demog aggregation"));
    await session.step(52, "And the table should have 4 rows", () => rowCount(page, 4));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
});
