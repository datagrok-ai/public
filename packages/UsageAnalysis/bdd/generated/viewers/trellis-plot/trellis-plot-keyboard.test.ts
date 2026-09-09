/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-keyboard.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cellsWideTall} from '../../../bindings/trellis-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKeyIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, noneOfFiltered} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, eventFired, listenFor, noErrors, propertyShouldBe, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot keyboard navigation", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-keyboard.feature", import.meta.url);
  test("Trellis plot keyboard navigation", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]));
    await session.step(16, "Then the cells of trellis plot viewer should be 2 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 2, 4));
    await run.scenario("Four arrows walk the current cell around a square and back", async () => {
      await session.step(19, "When user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(20, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(21, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(22, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(23, "When user presses ArrowDown in trellis plot viewer", () => pressKeyIn(page, "ArrowDown", el("trellis plot viewer")));
      await session.step(24, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Other\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Other"));
      await session.step(25, "When user presses ArrowLeft in trellis plot viewer", () => pressKeyIn(page, "ArrowLeft", el("trellis plot viewer")));
      await session.step(26, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Other\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Other"));
      await session.step(27, "When user presses ArrowUp in trellis plot viewer", () => pressKeyIn(page, "ArrowUp", el("trellis plot viewer")));
      await session.step(28, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Caucasian"));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every arrow announces the cell it moved to", async () => {
      await session.step(32, "Given user listens for \"d4-trellis-plot-current-cell-changed\" event on trellis plot viewer", () => listenFor(page, "d4-trellis-plot-current-cell-changed", el("trellis plot viewer")));
      await session.step(33, "When user clicks on the \"cell F | Asian\" area of trellis plot viewer", () => clickArea(page, "cell F | Asian", el("trellis plot viewer")));
      await session.step(34, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Asian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Asian"));
      await session.step(35, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(36, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Asian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Asian"));
      await session.step(37, "When user presses ArrowDown in trellis plot viewer", () => pressKeyIn(page, "ArrowDown", el("trellis plot viewer")));
      await session.step(38, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Black\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Black"));
      await session.step(39, "When user presses ArrowLeft in trellis plot viewer", () => pressKeyIn(page, "ArrowLeft", el("trellis plot viewer")));
      await session.step(40, "Then the \"current cell\" reading of trellis plot viewer should be \"F | Black\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "F | Black"));
      await session.step(41, "And \"d4-trellis-plot-current-cell-changed\" event should have fired on trellis plot viewer", () => eventFired(page, "d4-trellis-plot-current-cell-changed", el("trellis plot viewer")));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An arrow carries the filter to the cell it lands on", async () => {
      await session.step(45, "When user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(46, "And user clicks on the \"cell F | Caucasian\" area of trellis plot viewer", () => clickArea(page, "cell F | Caucasian", el("trellis plot viewer")));
      await session.step(47, "Then 480 rows should pass the filter", () => filterPasses(page, 480));
      await session.step(48, "When user presses ArrowRight in trellis plot viewer", () => pressKeyIn(page, "ArrowRight", el("trellis plot viewer")));
      await session.step(49, "Then the \"current cell\" reading of trellis plot viewer should be \"M | Caucasian\"", () => readingReads(page, "current cell", el("trellis plot viewer"), "M | Caucasian"));
      await session.step(50, "And 416 rows should pass the filter", () => filterPasses(page, 416));
      await session.step(51, "And no rows where \"SEX\" is \"F\" should pass the filter", () => noneOfFiltered(page, "SEX", "F"));
      await session.step(52, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Escape gives every row back", async () => {
      await session.step(55, "When user presses Escape in trellis plot viewer", () => pressKeyIn(page, "Escape", el("trellis plot viewer")));
      await session.step(56, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(57, "And the \"current cell\" reading of trellis plot viewer should be \"\"", () => readingReads(page, "current cell", el("trellis plot viewer"), ""));
      await session.step(58, "When user sets \"On Click\" property of trellis plot viewer to \"None\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(59, "Then \"On Click\" property of trellis plot viewer should be \"None\"", () => propertyShouldBe(page, "On Click", el("trellis plot viewer"), "None"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
