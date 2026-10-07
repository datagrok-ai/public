/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-viewer-filter-and-menus.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, pressKey, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, menuLists, noBalloons, noErrors, pickFromContextMenu, readingIs, readingNotAsRemembered, rememberReading, rightClickArea, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot viewer filter, menus and a second undo cycle", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-viewer-filter-and-menus.feature", import.meta.url);
  test("Trellis plot viewer filter, menus and a second undo cycle", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(18, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(19, "And the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
    await run.scenario("The viewer's Filter formula narrows its rows and leaves the table's filter alone", async () => {
      await session.step(22, "When user sets \"Filter\" property of trellis plot viewer to \"${AGE} > 40\"", () => setProperty(page, "Filter", el("trellis plot viewer"), "${AGE} > 40"));
      await session.step(23, "Then the \"rows shown\" reading of trellis plot viewer should be 635", () => readingIs(page, "rows shown", el("trellis plot viewer"), 635));
      await session.step(24, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(25, "When user sets \"Filter\" property of trellis plot viewer to \"\"", () => setProperty(page, "Filter", el("trellis plot viewer"), ""));
      await session.step(26, "Then the \"rows shown\" reading of trellis plot viewer should be 1000", () => readingIs(page, "rows shown", el("trellis plot viewer"), 1000));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A cell's context menu holds the inner viewer's group", async () => {
      await session.step(30, "When user right-clicks on the \"cell body F | Caucasian\" area of trellis plot viewer", () => rightClickArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(31, "Then the open menu should list \"Scatter plot > Lasso Tool\"", () => menuLists(page, "Scatter plot > Lasso Tool"));
      await session.step(32, "And the open menu should list \"Scatter plot > Markers\"", () => menuLists(page, "Scatter plot > Markers"));
      await session.step(33, "And the open menu should list \"Scatter plot > Selection\"", () => menuLists(page, "Scatter plot > Selection"));
      await session.step(34, "And the open menu should list \"General > Clone\"", () => menuLists(page, "General > Clone"));
      await session.step(35, "And the open menu should list \"Properties...\"", () => menuLists(page, "Properties..."));
      await session.step(36, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(37, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("To Script > To JavaScript prints the call that rebuilds the trellis", async () => {
      await session.step(40, "Given the package autostarts have completed", () => autostartsCompleted(page));
      await session.step(41, "When user picks \"To Script > To JavaScript\" from the context menu of trellis plot viewer", () => pickFromContextMenu(page, "To Script > To JavaScript", el("trellis plot viewer")));
      await session.step(42, "Then balloon should contain text \"addViewer\"", () => shouldContainText(page, el("balloon"), "addViewer"));
      await session.step(43, "And balloon should contain text \"Trellis\"", () => shouldContainText(page, el("balloon"), "Trellis"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The inner viewer's tab of the context panel changes every cell", async () => {
      await session.step(47, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(48, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(49, "When user clicks on \"Scatter plot\" tab in context panel", () => clickOn(page, el("\"Scatter plot\" tab in context panel")));
      await session.step(50, "Given \"X Axis\" category in context panel is expanded", () => isExpanded(page, el("\"X Axis\" category in context panel")));
      await session.step(51, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(52, "And user remembers the \"cell signature M | Asian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(53, "And user selects \"WEIGHT\" in \"X\" property in context panel", () => selectIn(page, "WEIGHT", el("\"X\" property in context panel")));
      await session.step(54, "Then \"xColumnName\" inner property of trellis plot viewer should be \"WEIGHT\"", () => innerPropertyShouldBe(page, "xColumnName", el("trellis plot viewer"), "WEIGHT"));
      await session.step(55, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(56, "And the \"cell signature M | Asian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Two cycles of undo and redo after the title-bar close leave no error", async () => {
      await session.step(60, "When user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
      await session.step(61, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(62, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(63, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(64, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(65, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(66, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(67, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(68, "Then the open tableview should have 1 trellis plot viewer", () => viewerCount(page, 1, "trellis plot"));
      await session.step(69, "And the \"cells drawn\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 8));
      await session.step(70, "When user presses Control+Shift+Z", () => pressKey(page, "Control+Shift+Z"));
      await session.step(71, "Then the open tableview should have 0 trellis plot viewers", () => viewerCount(page, 0, "trellis plot"));
      await session.step(72, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
