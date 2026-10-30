/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/table-properties.feature
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
import {clickOn, enterInto, expand, pressKey, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, dialogCloses, openTableOf, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A table made current shows its properties in the context panel, and its Rename... renames it", () => {
  const session = feature(test, "features/table-view/table-properties.feature", import.meta.url);
  test("A table made current shows its properties in the context panel, and its Rename... renames it", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(17, "And user opens a table \"tp\" with:", () => openTableOf(page, "tp", [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cy","Oslo","3"]]), [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cy","Oslo","3"]]);
    await session.step(22, "When user presses Alt+T", () => pressKey(page, "Alt+T"));
    await session.step(23, "And user clicks on the \"cell 1 of name\" area of Grid viewer in \"Tables\" dock panel", () => clickArea(page, "cell 1 of name", el("Grid viewer in \"Tables\" dock panel")));
    await session.step(24, "Then the context panel should show \"tp\"", () => contextPanelShows(page, "tp"));
    await run.scenario("General shows the table's row and column counts, and Actions its commands", async () => {
      await session.step(27, "When user expands General pane in context panel", () => expand(page, el("General pane in context panel")));
      await session.step(28, "Then General pane in context panel should contain the text \"Rows3\"", () => shouldContainText(page, el("General pane in context panel"), "Rows3"));
      await session.step(29, "And General pane in context panel should contain the text \"Columns3\"", () => shouldContainText(page, el("General pane in context panel"), "Columns3"));
      await session.step(30, "When user expands Actions pane in context panel", () => expand(page, el("Actions pane in context panel")));
      await session.step(31, "Then \"Rename...\" label in Actions pane in context panel should be visible", () => shouldBe(page, el("\"Rename...\" label in Actions pane in context panel"), "visible"));
      await session.step(32, "And \"Clone\" label in Actions pane in context panel should be visible", () => shouldBe(page, el("\"Clone\" label in Actions pane in context panel"), "visible"));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rename... in Actions renames the table, its view and its row in the Table Manager", async () => {
      await session.step(36, "When user expands Actions pane in context panel", () => expand(page, el("Actions pane in context panel")));
      await session.step(37, "And user clicks on \"Rename...\" label in Actions pane in context panel", () => clickOn(page, el("\"Rename...\" label in Actions pane in context panel")));
      await session.step(38, "Then \"Rename table\" dialog should be visible", () => shouldBe(page, el("\"Rename table\" dialog"), "visible"));
      await session.step(39, "And \"New name:\" input in \"Rename table\" dialog should have value \"tp\"", () => shouldHaveValue(page, el("\"New name:\" input in \"Rename table\" dialog"), "tp"));
      await session.step(40, "When user enters \"tp2\" into \"New name:\" input in \"Rename table\" dialog", () => enterInto(page, "tp2", el("\"New name:\" input in \"Rename table\" dialog")));
      await session.step(41, "And user clicks on OK button in \"Rename table\" dialog", () => clickOn(page, el("OK button in \"Rename table\" dialog")));
      await session.step(42, "Then the \"Rename table\" dialog should close", () => dialogCloses(page, "Rename table"));
      await session.step(43, "And table \"tp2\" should be open", () => tableOpen(page, "tp2"));
      await session.step(44, "And the \"tp2\" view should be current", () => viewIsCurrent(page, "tp2"));
      await session.step(45, "And the \"text of cell 1 of name\" reading of Grid viewer in \"Tables\" dock panel should be \"tp2\"", () => readingReads(page, "text of cell 1 of name", el("Grid viewer in \"Tables\" dock panel"), "tp2"));
      await session.step(46, "When user presses Alt+T", () => pressKey(page, "Alt+T"));
      await session.step(47, "Then \"Tables\" dock panel should be absent", () => shouldBe(page, el("\"Tables\" dock panel"), "absent"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
