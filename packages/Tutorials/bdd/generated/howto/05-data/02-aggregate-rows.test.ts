/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/02-aggregate-rows.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addToPivotRow} from '../../../bindings/howto-data.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clearSavedParameters} from '@datagrok-libraries/bdd/bindings/tiers/viewers/pivot-table';
import {clickArea, closeContextMenu, hoverArea, noErrors, pickFromAreaContextMenu, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/02-aggregate-rows.feature", import.meta.url);
  test("Aggregate rows", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "When user picks \"Data > Aggregate Rows...\" from the top menu", () => pickFromTopMenu(page, "Data > Aggregate Rows..."), undefined, "Aggregate Rows: summarize the table group by group");
    await session.step(22, "Then the open tableview should have 1 pivot table viewer", () => viewerCount(page, 1, "pivot table"));
    await session.step(23, "And the \"group by\" reading of pivot table viewer should be \"DIS_POP\"", () => readingReads(page, "group by", el("pivot table viewer"), "DIS_POP"));
    await session.step(25, "When user clicks on the \"remove group by chip DIS_POP\" area of pivot table viewer", () => clickArea(page, "remove group by chip DIS_POP", el("pivot table viewer")), undefined, "Group by: one result row per value");
    await session.step(26, "And user picks \"RACE\" with the plus of the \"group by\" row of pivot table viewer", () => addToPivotRow(page, "RACE", "group by"));
    await session.step(27, "Then the \"group by\" reading of pivot table viewer should be \"RACE\"", () => readingReads(page, "group by", el("pivot table viewer"), "RACE"));
    await session.step(29, "When user clicks on the \"remove pivot chip SEVERITY\" area of pivot table viewer", () => clickArea(page, "remove pivot chip SEVERITY", el("pivot table viewer")), undefined, "Pivot: one result column per value");
    await session.step(30, "And user picks \"SEX\" with the plus of the \"pivot\" row of pivot table viewer", () => addToPivotRow(page, "SEX", "pivot"));
    await session.step(31, "Then the \"pivot\" reading of pivot table viewer should be \"SEX\"", () => readingReads(page, "pivot", el("pivot table viewer"), "SEX"));
    await session.step(33, "When user picks \"WEIGHT\" with the plus of the \"aggregate\" row of pivot table viewer", () => addToPivotRow(page, "WEIGHT", "aggregate"), undefined, "Aggregate: what to compute for each group");
    await session.step(34, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), avg(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), avg(WEIGHT)"));
    await session.step(36, "When user picks \"Aggregation > med\" from the context menu of the \"aggregate chip avg(WEIGHT)\" area of pivot table viewer", () => pickFromAreaContextMenu(page, "Aggregation > med", "aggregate chip avg(WEIGHT)", el("pivot table viewer")), undefined, "Right-click an aggregate to change its function");
    await session.step(37, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(38, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(AGE), med(WEIGHT)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(AGE), med(WEIGHT)"));
    await session.step(39, "And the \"aggregated rows\" reading of pivot table viewer should be 4", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 4));
    await session.step(40, "When the video says \"ADD turns the result into a table of its own\" at the bottom", () => videoSaysAt(page, "ADD turns the result into a table of its own", "bottom"));
    await session.step(41, "And user hovers over the \"add to workspace\" area of pivot table viewer", () => hoverArea(page, "add to workspace", el("pivot table viewer")));
    await session.step(42, "And user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
    await session.step(43, "Then table \"demog aggregation\" should be open", () => tableOpen(page, "demog aggregation"));
    await session.step(44, "And the \"demog aggregation\" view should be current", () => viewIsCurrent(page, "demog aggregation"));
    await session.step(45, "And the table should have 4 rows", () => rowCount(page, 4));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
});
