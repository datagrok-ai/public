/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/08-end-to-end/04-clean-and-reshape.feature
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
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addToPivotRow} from '../../../bindings/howto-data.js';
import {pickAggregation} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, replaceCode, shouldBe, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRow, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clearSavedParameters} from '@datagrok-libraries/bdd/bindings/tiers/viewers/pivot-table';
import {clickArea, hoverArea, noErrors, painted, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · End to end", () => {
  const session = feature(test, "features/howto/08-end-to-end/04-clean-and-reshape.feature", import.meta.url);
  test("Clean and reshape a table", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await session.step(18, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "When the video says \"The task: the average BMI of each race, as a table and a chart\"", () => videoSays(page, "The task: the average BMI of each race, as a table and a chart"));
    await session.step(24, "And user picks \"Select > Missing Values...\" from the top menu", () => pickFromTopMenu(page, "Select > Missing Values..."), undefined, "Clean first: Select > Missing Values finds the rows with no height");
    await session.step(25, "Then Select Missing Values dialog should be visible", () => shouldBe(page, el("Select Missing Values dialog"), "visible"));
    await session.step(26, "When user toggles the \"HEIGHT\" column in the column list of Select Missing Values dialog", () => toggleInColumnList(page, "HEIGHT", el("Select Missing Values dialog")));
    await session.step(27, "And user clicks on OK button in Select Missing Values dialog", () => clickOn(page, el("OK button in Select Missing Values dialog")));
    await session.step(28, "Then 751 rows should be selected", () => selectedRowCount(page, 751));
    await session.step(30, "When user picks \"Edit > Remove > Selected Rows\" from the top menu", () => pickFromTopMenu(page, "Edit > Remove > Selected Rows"), undefined, "Edit > Remove > Selected Rows drops them");
    await session.step(31, "Then the table should have 5099 rows", () => rowCount(page, 5099));
    await session.step(33, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")), undefined, "Add New Column: BMI from weight and height");
    await session.step(34, "Then \"Add New Column\" dialog should be visible", () => shouldBe(page, el("\"Add New Column\" dialog"), "visible"));
    await session.step(35, "When user types \"BMI\" into column name input", () => typeInto(page, "BMI", el("column name input")));
    await session.step(36, "And user replaces the code of code editor in \"Add New Column\" dialog with \"${WEIGHT} / Pow(${HEIGHT} / 100, 2)\"", () => replaceCode(page, el("code editor in \"Add New Column\" dialog"), "${WEIGHT} / Pow(${HEIGHT} / 100, 2)"));
    await session.step(37, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
    await session.step(38, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
    await session.step(39, "And the table should have a column \"BMI\"", () => hasColumn(page, "BMI"));
    await session.step(41, "And the \"BMI\" cell of row 1 should be displayed as \"28.42\"", () => displayedInRow(page, "BMI", 1, "28.42"));
    await session.step(43, "When user picks \"Data > Aggregate Rows...\" from the top menu", () => pickFromTopMenu(page, "Data > Aggregate Rows..."), undefined, "Reshape: Data > Aggregate Rows summarizes group by group");
    await session.step(44, "Then the open tableview should have 1 pivot table viewer", () => viewerCount(page, 1, "pivot table"));
    await session.step(46, "When user clicks on the \"remove group by chip DIS_POP\" area of pivot table viewer", () => clickArea(page, "remove group by chip DIS_POP", el("pivot table viewer")), undefined, "One row per race");
    await session.step(47, "And user picks \"RACE\" with the plus of the \"group by\" row of pivot table viewer", () => addToPivotRow(page, "RACE", "group by"));
    await session.step(48, "And user clicks on the \"remove pivot chip SEVERITY\" area of pivot table viewer", () => clickArea(page, "remove pivot chip SEVERITY", el("pivot table viewer")));
    await session.step(49, "Then the \"group by\" reading of pivot table viewer should be \"RACE\"", () => readingReads(page, "group by", el("pivot table viewer"), "RACE"));
    await session.step(51, "When user clicks on the \"remove aggregate chip avg(AGE)\" area of pivot table viewer", () => clickArea(page, "remove aggregate chip avg(AGE)", el("pivot table viewer")), undefined, "with the average BMI of each");
    await session.step(52, "And user picks \"BMI\" with the plus of the \"aggregate\" row of pivot table viewer", () => addToPivotRow(page, "BMI", "aggregate"));
    await session.step(53, "Then the \"aggregate\" reading of pivot table viewer should be \"avg(BMI)\"", () => readingReads(page, "aggregate", el("pivot table viewer"), "avg(BMI)"));
    await session.step(54, "And the \"aggregated rows\" reading of pivot table viewer should be 4", () => readingIs(page, "aggregated rows", el("pivot table viewer"), 4));
    await session.step(56, "When user hovers over the \"add to workspace\" area of pivot table viewer", () => hoverArea(page, "add to workspace", el("pivot table viewer")), undefined, "ADD makes the summary a table of its own");
    await session.step(57, "And user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
    await session.step(58, "Then the \"demog aggregation\" view should be current", () => viewIsCurrent(page, "demog aggregation"));
    await session.step(59, "And the table should have 4 rows", () => rowCount(page, 4));
    await session.step(61, "When user clicks on \"bar chart\" icon in toolbox", () => clickOn(page, el("\"bar chart\" icon in toolbox")), undefined, "Chart it: a bar chart of the average BMI per race");
    await session.step(62, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(63, "When user picks \"avg(BMI)\" in the \"value\" column selector of bar chart viewer", () => pickInColumnSelector(page, "avg(BMI)", "value", el("bar chart viewer")));
    await session.step(64, "And user picks the \"avg\" aggregation in the \"value\" column selector of bar chart viewer", () => pickAggregation(page, "avg", "value", el("bar chart viewer")));
    await session.step(65, "Then bar chart viewer should be painted", () => painted(page, el("bar chart viewer")));
    await session.step(66, "When the video says \"From raw rows to an answer: four numbers and a chart\"", () => videoSays(page, "From raw rows to an answer: four numbers and a chart"));
    await session.step(67, "And user hovers over bar chart viewer", () => hoverOver(page, el("bar chart viewer")));
    await session.step(68, "Then no errors should have been logged", () => noErrors(page));
  });
});
