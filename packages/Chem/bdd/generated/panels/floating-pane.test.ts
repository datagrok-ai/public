/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/floating-pane.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {expand, hoverText, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, closePanelWatched, closedPanelReleased, contextPanelOpen, contextPanelShows, contextPanelShowsCurrentCell, dragPaneOut, openDataset, openDatasetRowsAs, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, loadLayout, noErrors, saveLayout} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dockToViewEdge} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A pane dragged out of the context panel becomes a viewer of the current molecule", () => {
  const session = feature(test, "features/panels/floating-pane.feature", import.meta.url);
  test("A pane dragged out of the context panel becomes a viewer of the current molecule", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(20, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await session.step(21, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(22, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
    await session.step(23, "Then row 1 should be current", () => currentRowIs(page, 1));
    await session.step(24, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
    await session.step(25, "When user expands Chemistry accordion header in context panel", () => expand(page, el("Chemistry accordion header in context panel")));
    await session.step(26, "And user drags Properties accordion header in context panel out of the accordion", () => dragPaneOut(page, el("Properties accordion header in context panel")));
    await session.step(27, "Then \"Properties\" floating panel should contain the text \"C24H23NO3\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C24H23NO3"));
    await run.scenario("The dragged-out pane recalculates for each molecule", async () => {
      await session.step(30, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
      await session.step(31, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(32, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(33, "And \"Properties\" floating panel should contain the text \"C12H14\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C12H14"));
      await session.step(34, "And \"Properties\" floating panel should not contain the text \"C24H23NO3\"", () => shouldNotContainText(page, el("\"Properties\" floating panel"), "C24H23NO3"));
      await session.step(35, "When user clicks on the \"cell 6 of canonical_smiles\" area of grid", () => clickArea(page, "cell 6 of canonical_smiles", el("grid")));
      await session.step(36, "Then row 6 should be current", () => currentRowIs(page, 6));
      await session.step(37, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(38, "And \"Properties\" floating panel should contain the text \"C27H35NO3\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C27H35NO3"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Anything else made current has no effect on it", async () => {
      await session.step(42, "When user clicks on the \"cell 6 of canonical_smiles\" area of grid", () => clickArea(page, "cell 6 of canonical_smiles", el("grid")));
      await session.step(43, "Then row 6 should be current", () => currentRowIs(page, 6));
      await session.step(44, "And \"Properties\" floating panel should contain the text \"C27H35NO3\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C27H35NO3"));
      await session.step(45, "When user clicks on the \"cell 3 of molregno\" area of grid", () => clickArea(page, "cell 3 of molregno", el("grid")));
      await session.step(46, "Then row 3 should be current", () => currentRowIs(page, 3));
      await session.step(47, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(48, "And \"Properties\" floating panel should contain the text \"C27H35NO3\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C27H35NO3"));
      await session.step(49, "When user clicks on the \"header molregno\" area of grid", () => clickArea(page, "header molregno", el("grid")));
      await session.step(50, "Then the context panel should show \"molregno\"", () => contextPanelShows(page, "molregno"));
      await session.step(51, "And \"Properties\" floating panel should contain the text \"C27H35NO3\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C27H35NO3"));
      await session.step(52, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(53, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(54, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(55, "And \"Properties\" floating panel should contain the text \"C24H23NO3\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C24H23NO3"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("It comes back with the layout while it floats", async () => {
      await session.step(59, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
      await session.step(60, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(61, "And \"Properties\" floating panel should contain the text \"C12H14\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C12H14"));
      await session.step(62, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(63, "And user loads the saved layout", () => loadLayout(page));
      await session.step(64, "Then \"Properties\" floating panel should contain the text \"C12H14\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C12H14"));
      await session.step(65, "When user clicks on the \"cell 6 of canonical_smiles\" area of grid", () => clickArea(page, "cell 6 of canonical_smiles", el("grid")));
      await session.step(66, "Then row 6 should be current", () => currentRowIs(page, 6));
      await session.step(67, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(68, "And \"Properties\" floating panel should contain the text \"C27H35NO3\"", () => shouldContainText(page, el("\"Properties\" floating panel"), "C27H35NO3"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("It docks into the view, follows there, and comes back with the layout docked", async () => {
      await session.step(72, "When user docks info panel viewer to the right edge of the view", () => dockToViewEdge(page, el("info panel viewer"), "right"));
      await session.step(73, "Then \"Properties\" floating panel should be absent", () => shouldBe(page, el("\"Properties\" floating panel"), "absent"));
      await session.step(74, "And \"Properties\" dock panel should contain the text \"C27H35NO3\"", () => shouldContainText(page, el("\"Properties\" dock panel"), "C27H35NO3"));
      await session.step(75, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
      await session.step(76, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(77, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(78, "And \"Properties\" dock panel should contain the text \"C12H14\"", () => shouldContainText(page, el("\"Properties\" dock panel"), "C12H14"));
      await session.step(79, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(80, "And user loads the saved layout", () => loadLayout(page));
      await session.step(81, "Then \"Properties\" floating panel should be absent", () => shouldBe(page, el("\"Properties\" floating panel"), "absent"));
      await session.step(82, "And \"Properties\" dock panel should contain the text \"C12H14\"", () => shouldContainText(page, el("\"Properties\" dock panel"), "C12H14"));
      await session.step(83, "When user clicks on the \"cell 1 of canonical_smiles\" area of grid", () => clickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(84, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(85, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(86, "And \"Properties\" dock panel should contain the text \"C24H23NO3\"", () => shouldContainText(page, el("\"Properties\" dock panel"), "C24H23NO3"));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Out of sight it waits, and shows the latest molecule once it is back", async () => {
      await session.step(90, "Given user opens smiles-50 dataset keeping the first 10 rows as \"more molecules\"", () => openDatasetRowsAs(page, ds("smiles-50"), 10, "more molecules"));
      await session.step(91, "And the context panel is open", () => contextPanelOpen(page));
      await session.step(92, "When user clicks on the \"cell 6 of canonical_smiles\" area of grid", () => clickArea(page, "cell 6 of canonical_smiles", el("grid")));
      await session.step(93, "Then row 6 should be current", () => currentRowIs(page, 6));
      await session.step(94, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(95, "When user switches to the \"smiles-50\" table view", () => switchTableView(page, "smiles-50"));
      await session.step(96, "Then \"Properties\" dock panel should contain the text \"C27H35NO3\"", () => shouldContainText(page, el("\"Properties\" dock panel"), "C27H35NO3"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Closing it ends it, and nothing keeps it", async () => {
      await session.step(100, "When user hovers over the text \"Molecular formula\" in \"Properties\" dock panel", () => hoverText(page, "Molecular formula", el("\"Properties\" dock panel")));
      await session.step(101, "And user closes \"Properties\" dock panel, watching that nothing keeps it", () => closePanelWatched(page, el("\"Properties\" dock panel")));
      await session.step(102, "Then info panel viewer should be absent", () => shouldBe(page, el("info panel viewer"), "absent"));
      await session.step(103, "When user clicks on the \"cell 5 of canonical_smiles\" area of grid", () => clickArea(page, "cell 5 of canonical_smiles", el("grid")));
      await session.step(104, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(105, "And the context panel should show the current cell", () => contextPanelShowsCurrentCell(page));
      await session.step(106, "And info panel viewer should be absent", () => shouldBe(page, el("info panel viewer"), "absent"));
      await session.step(107, "And nothing should keep the closed panel", () => closedPanelReleased(page));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
