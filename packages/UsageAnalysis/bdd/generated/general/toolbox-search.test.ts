/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/toolbox-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The table view's Search box filters the rows its text matches", () => {
  const session = feature(test, "features/general/toolbox-search.feature", import.meta.url);
  test("The table view's Search box filters the rows its text matches", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(22, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(23, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(24, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await run.scenario("A condition on a numeric column filters the rows it names", async () => {
      await session.step(27, "When user types \"AGE > 50\" into table search", () => typeInto(page, "AGE > 50", el("table search")));
      await session.step(28, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(29, "Then 2176 rows should pass the filter", () => filterPasses(page, 2176));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A condition on a text column filters the rows it names", async () => {
      await session.step(33, "When user types \"SEX = M\" into table search", () => typeInto(page, "SEX = M", el("table search")));
      await session.step(34, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(35, "Then 2607 rows should pass the filter", () => filterPasses(page, 2607));
      await session.step(36, "When user types \"RACE = Asian\" into table search", () => typeInto(page, "RACE = Asian", el("table search")));
      await session.step(37, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(38, "Then 72 rows should pass the filter", () => filterPasses(page, 72));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An empty Search box gives every row back", async () => {
      await session.step(42, "When user clears table search", () => clearField(page, el("table search")));
      await session.step(43, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(44, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
