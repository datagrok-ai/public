/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/table-manager.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, pressKey, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, contextPanelOpen, contextPanelShows, dialogCloses, openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, noErrors, pickFromAreaContextMenu, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Table Manager lists the open tables and moves between them", () => {
  const session = feature(test, "features/general/table-manager.feature", import.meta.url);
  test("The Table Manager lists the open tables and moves between them", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(19, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(20, "And user opens iris dataset", () => openDataset(page, ds("iris")));
    await session.step(21, "And user opens beer dataset", () => openDataset(page, ds("beer")));
    await session.step(22, "Then the open table views should be exactly \"cars, iris, beer\"", () => openTableViewsExactly(page, "cars, iris, beer"));
    await run.scenario("Alt+T docks the manager with every open table, in the order they opened", async () => {
      await session.step(25, "Then \"Tables\" dock panel should be absent", () => shouldBe(page, el("\"Tables\" dock panel"), "absent"));
      await session.step(26, "When user presses Alt+T", () => pressKey(page, "Alt+T"));
      await session.step(27, "Then the \"rows\" reading of Grid viewer in \"Tables\" dock panel should be 3", () => readingIs(page, "rows", el("Grid viewer in \"Tables\" dock panel"), 3));
      await session.step(28, "And the \"text of cell 1 of name\" reading of Grid viewer in \"Tables\" dock panel should be \"cars\"", () => readingReads(page, "text of cell 1 of name", el("Grid viewer in \"Tables\" dock panel"), "cars"));
      await session.step(29, "And the \"text of cell 2 of name\" reading of Grid viewer in \"Tables\" dock panel should be \"iris\"", () => readingReads(page, "text of cell 2 of name", el("Grid viewer in \"Tables\" dock panel"), "iris"));
      await session.step(30, "And the \"text of cell 3 of name\" reading of Grid viewer in \"Tables\" dock panel should be \"beer\"", () => readingReads(page, "text of cell 3 of name", el("Grid viewer in \"Tables\" dock panel"), "beer"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a row brings its table's view to front and makes the table current", async () => {
      await session.step(34, "When user clicks on the \"cell 1 of name\" area of Grid viewer in \"Tables\" dock panel", () => clickArea(page, "cell 1 of name", el("Grid viewer in \"Tables\" dock panel")));
      await session.step(35, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
      await session.step(36, "And the context panel should show \"cars\"", () => contextPanelShows(page, "cars"));
      await session.step(37, "When user clicks on the \"cell 2 of name\" area of Grid viewer in \"Tables\" dock panel", () => clickArea(page, "cell 2 of name", el("Grid viewer in \"Tables\" dock panel")));
      await session.step(38, "Then the \"iris\" view should be current", () => viewIsCurrent(page, "iris"));
      await session.step(39, "And the context panel should show \"iris\"", () => contextPanelShows(page, "iris"));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rename... in the current table's Actions pane renames the table, its view and its row in the manager", async () => {
      await session.step(43, "When user expands Actions pane in context panel", () => expand(page, el("Actions pane in context panel")));
      await session.step(44, "And user clicks on \"Rename...\" label in Actions pane in context panel", () => clickOn(page, el("\"Rename...\" label in Actions pane in context panel")));
      await session.step(45, "Then \"Rename table\" dialog should be visible", () => shouldBe(page, el("\"Rename table\" dialog"), "visible"));
      await session.step(46, "And \"New name:\" input in \"Rename table\" dialog should have value \"iris\"", () => shouldHaveValue(page, el("\"New name:\" input in \"Rename table\" dialog"), "iris"));
      await session.step(47, "When user enters \"iris2\" into \"New name:\" input in \"Rename table\" dialog", () => enterInto(page, "iris2", el("\"New name:\" input in \"Rename table\" dialog")));
      await session.step(48, "And user clicks on OK button in \"Rename table\" dialog", () => clickOn(page, el("OK button in \"Rename table\" dialog")));
      await session.step(49, "Then the \"Rename table\" dialog should close", () => dialogCloses(page, "Rename table"));
      await session.step(50, "And the \"iris2\" view should be current", () => viewIsCurrent(page, "iris2"));
      await session.step(51, "And the \"text of cell 2 of name\" reading of Grid viewer in \"Tables\" dock panel should be \"iris2\"", () => readingReads(page, "text of cell 2 of name", el("Grid viewer in \"Tables\" dock panel"), "iris2"));
      await session.step(52, "When user clicks on \"Rename...\" label in Actions pane in context panel", () => clickOn(page, el("\"Rename...\" label in Actions pane in context panel")));
      await session.step(53, "And user enters \"iris\" into \"New name:\" input in \"Rename table\" dialog", () => enterInto(page, "iris", el("\"New name:\" input in \"Rename table\" dialog")));
      await session.step(54, "And user clicks on OK button in \"Rename table\" dialog", () => clickOn(page, el("OK button in \"Rename table\" dialog")));
      await session.step(55, "Then the \"Rename table\" dialog should close", () => dialogCloses(page, "Rename table"));
      await session.step(56, "And the \"text of cell 2 of name\" reading of Grid viewer in \"Tables\" dock panel should be \"iris\"", () => readingReads(page, "text of cell 2 of name", el("Grid viewer in \"Tables\" dock panel"), "iris"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Open as table makes a table of the manager's list, and a closed table leaves the list", async () => {
      await session.step(60, "When user picks \"Open as table\" from the context menu of the \"cell 1 of name\" area of Grid viewer in \"Tables\" dock panel", () => pickFromAreaContextMenu(page, "Open as table", "cell 1 of name", el("Grid viewer in \"Tables\" dock panel")));
      await session.step(61, "Then the \"rows\" reading of Grid viewer in \"Tables\" dock panel should be 4", () => readingIs(page, "rows", el("Grid viewer in \"Tables\" dock panel"), 4));
      await session.step(62, "And the table should have 3 rows", () => rowCount(page, 3));
      await session.step(63, "And the value of \"name\" column in row 1 should be \"cars\"", () => valueInRow(page, "name", 1, "cars"));
      await session.step(64, "And the value of \"name\" column in row 3 should be \"beer\"", () => valueInRow(page, "name", 3, "beer"));
      await session.step(65, "When user closes the current view", () => closeCurrentView(page));
      await session.step(66, "Then the open table views should be exactly \"cars, iris, beer\"", () => openTableViewsExactly(page, "cars, iris, beer"));
      await session.step(67, "And the \"rows\" reading of Grid viewer in \"Tables\" dock panel should be 3", () => readingIs(page, "rows", el("Grid viewer in \"Tables\" dock panel"), 3));
      await session.step(68, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Alt+T closes the manager again", async () => {
      await session.step(71, "When user presses Alt+T", () => pressKey(page, "Alt+T"));
      await session.step(72, "Then \"Tables\" dock panel should be absent", () => shouldBe(page, el("\"Tables\" dock panel"), "absent"));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
