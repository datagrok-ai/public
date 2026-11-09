/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/general/toolbox-search.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [public#3963]
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
import {clearField, clickOn, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterPasses, filterPassesAll, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The table view's Search box filters the rows its text matches", () => {
  const session = feature(test, "features/general/toolbox-search.feature", import.meta.url);
  test("The table view's Search box filters the rows its text matches", {tag: ["@journey", "@realizes:public#3963"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(24, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await run.scenario("A condition on a numeric column filters the rows it names", async () => {
      await session.step(27, "When user types \"AGE > 50\" into table search", () => typeInto(page, "AGE > 50", el("table search")));
      await session.step(28, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(29, "Then 367 rows should pass the filter", () => filterPasses(page, 367));
      await session.step(30, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A condition on a text column filters the rows it names", async () => {
      await session.step(33, "When user types \"SEX = M\" into table search", () => typeInto(page, "SEX = M", el("table search")));
      await session.step(34, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(35, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(36, "When user types \"RACE = Asian\" into table search", () => typeInto(page, "RACE = Asian", el("table search")));
      await session.step(37, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(38, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An empty Search box gives every row back, and Extract Rows takes only the rows a search left", async () => {
      await session.step(43, "When user clears table search", () => clearField(page, el("table search")));
      await session.step(44, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(45, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
      await session.step(47, "When user types \"RACE = Asian\" into table search", () => typeInto(page, "RACE = Asian", el("table search")));
      await session.step(48, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(49, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(50, "When user clicks on \"Filtered: 15\" text in status bar", () => clickOn(page, el("\"Filtered: 15\" text in status bar")));
      await session.step(51, "And user picks \"Extract Rows\" from the open menu", () => pickFromOpenMenu(page, "Extract Rows"));
      await session.step(52, "Then the table should have 15 rows", () => rowCount(page, 15));
      await session.step(53, "And every value of \"RACE\" column should match \"^Asian$\"", () => everyValueMatches(page, "RACE", "^Asian$"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
