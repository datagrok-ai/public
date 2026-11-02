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
import {pressKeyIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, noneOfFiltered} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, boundTable, clickArea, eventFired, listenFor, noErrors, propertyShouldBe, readingReads, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot keyboard navigation", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-keyboard.feature", import.meta.url);
  test("Trellis plot keyboard navigation", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot", "@realizes:core#1057"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(18, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await run.scenario("Four arrows walk the current cell around a square and back", async () => {
      await session.step(21, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(22, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(23, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(24, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(25, "When user presses ArrowDown in trellis plot viewer", () => pressKeyIn(page, "ArrowDown", el("trellis plot viewer")));
      await session.step(26, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Other\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Other"));
      await session.step(27, "When user presses ArrowLeft in trellis plot viewer", () => pressKeyIn(page, "ArrowLeft", el("trellis plot viewer")));
      await session.step(28, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Other\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Other"));
      await session.step(29, "When user presses ArrowUp in trellis plot viewer", () => pressKeyIn(page, "ArrowUp", el("trellis plot viewer")));
      await session.step(30, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every arrow announces the cell it moved to", async () => {
      await session.step(34, "Given user listens for \"d4-trellis-plot-current-cell-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-current-cell-changed", el("trellis plot viewer")));
      await session.step(35, "When user clicks on the \"cell F | Asian\" area of trellis plot viewer", () => clickArea(page, "cell F | Asian", el("trellis plot viewer")));
      await session.step(36, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Asian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Asian"));
      await session.step(37, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(38, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Asian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Asian"));
      await session.step(39, "When user presses ArrowDown in trellis plot viewer", () => pressKeyIn(page, "ArrowDown", el("trellis plot viewer")));
      await session.step(40, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Black\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Black"));
      await session.step(41, "When user presses ArrowLeft in trellis plot viewer", () => pressKeyIn(page, "ArrowLeft", el("trellis plot viewer")));
      await session.step(42, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Black\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Black"));
      await session.step(43, "And \"d4-trellis-plot-current-cell-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-current-cell-changed", el("trellis plot viewer")));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An arrow carries the filter to the cell it lands on", async () => {
      await session.step(47, "When user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(48, "And user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(49, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(50, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(51, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(52, "And 416 rows should pass the filter", () => filterPasses(page, 416));
      await session.step(53, "And no rows where \"SEX\" is \"F\" should pass the filter", () => noneOfFiltered(page, "SEX", "F"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape gives every row back", async () => {
      await session.step(58, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(59, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(60, "And the \"current cell\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "current cell", el("trellis plot viewer"), ""));
      await session.step(61, "When user sets \"On Click\" property of trellis plot viewer to \"None\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(62, "Then \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
      await session.step(64, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(65, "And user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(66, "When user sets \"table\" property of trellis plot viewer to \"cars\"", () => setProperty(page, "table", el("trellis plot viewer"), "cars"));
      await session.step(67, "Then trellis plot viewer should be bound to table \"cars\"", () => boundTable(page, el("trellis plot viewer"), "cars"));
      await session.step(68, "When user sets \"table\" property of trellis plot viewer to \"demog-1000\"", () => setProperty(page, "table", el("trellis plot viewer"), "demog-1000"));
      await session.step(69, "And user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
      await session.step(73, "Then trellis plot viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("trellis plot viewer"), "demog-1000"));
      await session.step(74, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(75, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(76, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(77, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
