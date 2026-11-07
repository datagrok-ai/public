/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-keyboard.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot, core#1057]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, noneOfFiltered} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, clickArea, eventFired, listenFor, noErrors, propertyShouldBe, readingReads, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot keyboard navigation", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-keyboard.feature", import.meta.url);
  test("Trellis plot keyboard navigation", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot", "@edge-cases", "@realizes:core#1057"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(21, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await run.scenario("Four arrows walk the current cell around a square and back", async () => {
      await session.step(24, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(25, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(26, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(27, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(28, "When user presses ArrowDown in trellis plot viewer", () => pressKeyIn(page, "ArrowDown", el("trellis plot viewer")));
      await session.step(29, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Other\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Other"));
      await session.step(30, "When user presses ArrowLeft in trellis plot viewer", () => pressKeyIn(page, "ArrowLeft", el("trellis plot viewer")));
      await session.step(31, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Other\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Other"));
      await session.step(32, "When user presses ArrowUp in trellis plot viewer", () => pressKeyIn(page, "ArrowUp", el("trellis plot viewer")));
      await session.step(33, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every arrow announces the cell it moved to", async () => {
      await session.step(37, "Given user listens for \"d4-trellis-plot-current-cell-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-current-cell-changed", el("trellis plot viewer")));
      await session.step(38, "When user clicks on the \"cell F | Asian\" area of trellis plot viewer", () => clickArea(page, "cell F | Asian", el("trellis plot viewer")));
      await session.step(39, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Asian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Asian"));
      await session.step(40, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(41, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Asian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Asian"));
      await session.step(42, "When user presses ArrowDown in trellis plot viewer", () => pressKeyIn(page, "ArrowDown", el("trellis plot viewer")));
      await session.step(43, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Black\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Black"));
      await session.step(44, "When user presses ArrowLeft in trellis plot viewer", () => pressKeyIn(page, "ArrowLeft", el("trellis plot viewer")));
      await session.step(45, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Black\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Black"));
      await session.step(46, "And \"d4-trellis-plot-current-cell-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-current-cell-changed", el("trellis plot viewer")));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An arrow carries the filter to the cell it lands on", async () => {
      await session.step(50, "When user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(51, "And user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(52, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(53, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(54, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(55, "And 416 rows should pass the filter", () => filterPasses(page, 416));
      await session.step(56, "And no rows where \"SEX\" is \"F\" should pass the filter", () => noneOfFiltered(page, "SEX", "F"));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape gives every row back", async () => {
      await session.step(60, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(61, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(62, "And the \"current cell\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "current cell", el("trellis plot viewer"), ""));
      await session.step(63, "When user sets \"On Click\" property of trellis plot viewer to \"None\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(64, "Then \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("After the trellis is set to another table and back, the arrow keys still move the current cell (core#1057)", async () => {
      await session.step(69, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(70, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(71, "And the context panel is open", () => contextPanelOpen(page));
      await session.step(72, "When user clicks on settings icon of trellis plot viewer", () => clickOn(page, el("settings icon of trellis plot viewer")));
      await session.step(73, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(74, "When user selects \"cars\" in \"Table\" property in context panel", () => selectIn(page, "cars", el("\"Table\" property in context panel")));
      await session.step(75, "Then trellis plot viewer should be bound to table \"cars\"", () => boundTable(page, el("trellis plot viewer"), "cars"));
      await session.step(76, "When user selects \"demog-1000\" in \"Table\" property in context panel", () => selectIn(page, "demog-1000", el("\"Table\" property in context panel")));
      await session.step(77, "Then trellis plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("trellis plot viewer"), "demog-1000"));
      await session.step(78, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
      await session.step(82, "And user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(83, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(84, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(85, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(86, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
