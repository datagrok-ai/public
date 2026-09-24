/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/numeric-matching.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '../../../bindings/grid.js';
import '../../../bindings/nx.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, pressKeyIn, selectIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactly, filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickPanelMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Numeric matching in the search box and the expression filter", () => {
  const session = feature(test, "features/viewers/filter-panel/numeric-matching.feature", import.meta.url);
  test("Numeric matching in the search box and the expression filter", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(15, "And user opens a table \"scores\" with:", () => openTableOf(page, "scores", [["Score","Prediction"],["0.47","Medium"],["0.5","Low"],["0.47","High"],["0.83","Medium"],["1.25","Low"]]), [["Score","Prediction"],["0.47","Medium"],["0.5","Low"],["0.47","High"],["0.83","Medium"],["1.25","Low"]]);
    await session.step(22, "Then all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("Ctrl+F, a float typed as the cell shows it, Enter — the search box filters the rows", async () => {
      await session.step(25, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
      await session.step(26, "Then table search should be visible", () => shouldBe(page, el("table search"), "visible"));
      await session.step(27, "When user types \"0.47\" into table search", () => typeInto(page, "0.47", el("table search")));
      await session.step(28, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(29, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
      await session.step(30, "When user types \">= 0.83\" into table search", () => typeInto(page, ">= 0.83", el("table search")));
      await session.step(31, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(32, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
      await session.step(33, "When user types \"< 0.47\" into table search", () => typeInto(page, "< 0.47", el("table search")));
      await session.step(34, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(35, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(36, "When user types \"0.47-0.5\" into table search", () => typeInto(page, "0.47-0.5", el("table search")));
      await session.step(37, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(38, "Then 3 rows should pass the filter", () => filterPasses(page, 3));
      await session.step(39, "When user clears table search", () => clearField(page, el("table search")));
      await session.step(40, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(41, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The expression filter compares a float column against the value it shows", async () => {
      await session.step(45, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(46, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(47, "And user selects \"Score\" in Column input in \"Expression\" filter card", () => selectIn(page, "Score", el("Column input in \"Expression\" filter card")));
      await session.step(48, "And user selects \"=\" in Operation input in \"Expression\" filter card", () => selectIn(page, "=", el("Operation input in \"Expression\" filter card")));
      await session.step(49, "And user types \"0.47\" into Value input in \"Expression\" filter card", () => typeInto(page, "0.47", el("Value input in \"Expression\" filter card")));
      await session.step(50, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(51, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
      await session.step(52, "And the \"categories of Expression\" reading of filter panel should be \"${Score} = 0.47\"", () => readingReads(page, "categories of Expression", el("filter panel"), "${Score} = 0.47"));
      await session.step(53, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(54, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(55, "And user selects \"Score\" in Column input in \"Expression\" filter card", () => selectIn(page, "Score", el("Column input in \"Expression\" filter card")));
      await session.step(56, "And user selects \">=\" in Operation input in \"Expression\" filter card", () => selectIn(page, ">=", el("Operation input in \"Expression\" filter card")));
      await session.step(57, "And user types \"0.83\" into Value input in \"Expression\" filter card", () => typeInto(page, "0.83", el("Value input in \"Expression\" filter card")));
      await session.step(58, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(59, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
      await session.step(60, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(61, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(62, "And user selects \"All Columns\" in Column input in \"Expression\" filter card", () => selectIn(page, "All Columns", el("Column input in \"Expression\" filter card")));
      await session.step(63, "And user selects \"equals\" in Operation input in \"Expression\" filter card", () => selectIn(page, "equals", el("Operation input in \"Expression\" filter card")));
      await session.step(64, "And user types \"0.47\" into Value input in \"Expression\" filter card", () => typeInto(page, "0.47", el("Value input in \"Expression\" filter card")));
      await session.step(65, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(66, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
      await session.step(67, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(68, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("All Columns with the word \"equals\" reaches the numeric columns of spgi", async () => {
      await session.step(72, "When user opens spgi dataset", () => openDataset(page, ds("spgi")));
      await session.step(73, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(74, "And user picks \"Add Filter | Expression\" from the filter panel menu", () => pickPanelMenu(page, "Add Filter | Expression"));
      await session.step(75, "And user selects \"All Columns\" in Column input in \"Expression\" filter card", () => selectIn(page, "All Columns", el("Column input in \"Expression\" filter card")));
      await session.step(76, "And user selects \"equals\" in Operation input in \"Expression\" filter card", () => selectIn(page, "equals", el("Operation input in \"Expression\" filter card")));
      await session.step(77, "And user types \"634783\" into Value input in \"Expression\" filter card", () => typeInto(page, "634783", el("Value input in \"Expression\" filter card")));
      await session.step(78, "And user clicks on \"Add filter\" button in \"Expression\" filter card", () => clickOn(page, el("\"Add filter\" button in \"Expression\" filter card")));
      await session.step(79, "Then 1 row should pass the filter", () => filterPasses(page, 1));
      await session.step(80, "And the filter should pass exactly the rows where \"CAST Idea ID\" is between 634783 and 634783", () => filterIsExactly(page, "CAST Idea ID", 634783, 634783));
      await session.step(81, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(82, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
