/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/file-drop-overlay.feature
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
import {pointerOver, pressKey, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, fileDragEndLost, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {clickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The file drop overlay always goes away", () => {
  const session = feature(test, "features/general/file-drop-overlay.feature", import.meta.url);
  test("The file drop overlay always goes away", {tag: ["@journey", "@GROK-20507"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await run.scenario("Escape takes the layer away while the grid has the focus", async () => {
      await session.step(21, "When user clicks on the \"cell 2 of model\" area of grid", () => clickArea(page, "cell 2 of model", el("grid")));
      await session.step(22, "And a file dragged from disk enters grid and the drag's end never reaches the page", () => fileDragEndLost(page, el("grid")));
      await session.step(23, "Then file drop overlay should be visible", () => shouldBe(page, el("file drop overlay"), "visible"));
      await session.step(24, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(25, "Then file drop overlay should be absent", () => shouldBe(page, el("file drop overlay"), "absent"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pointer move takes the layer away", async () => {
      await session.step(29, "When a file dragged from disk enters grid and the drag's end never reaches the page", () => fileDragEndLost(page, el("grid")));
      await session.step(30, "Then file drop overlay should be visible", () => shouldBe(page, el("file drop overlay"), "visible"));
      await session.step(31, "When user moves the pointer over file drop overlay", () => pointerOver(page, el("file drop overlay")));
      await session.step(32, "Then file drop overlay should be absent", () => shouldBe(page, el("file drop overlay"), "absent"));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A file dropped after the layer was dismissed still opens", async () => {
      await session.step(36, "When user drops the \"fixtures/cars-small.csv\" file of the project onto grid", () => dropFile(page, "fixtures/cars-small.csv", el("grid")));
      await session.step(37, "Then the \"cars-small\" table view should open with 3 rows", () => tableViewOpened(page, "cars-small", 3));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
      await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
