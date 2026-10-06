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
import {columnIncomplete, columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("The table view's Search box filters by compound and date conditions", () => {
  const session = feature(test, "features/general/toolbox-search.feature", import.meta.url);
  test("A single condition filters the rows it names", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(18, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await session.step(21, "When user types \"AGE > 50\" into table search", () => typeInto(page, "AGE > 50", el("table search")));
    await session.step(22, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(23, "Then 2176 rows should pass the filter", () => filterPasses(page, 2176));
    await session.step(24, "When user types \"SEX = M\" into table search", () => typeInto(page, "SEX = M", el("table search")));
    await session.step(25, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(26, "Then 2607 rows should pass the filter", () => filterPasses(page, 2607));
    await session.step(27, "When user types \"CONTROL = true\" into table search", () => typeInto(page, "CONTROL = true", el("table search")));
    await session.step(28, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(29, "Then 39 rows should pass the filter", () => filterPasses(page, 39));
    await session.step(31, "When user types \"STARTED > 1/1/1990\" into table search", () => typeInto(page, "STARTED > 1/1/1990", el("table search")));
    await session.step(32, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(33, "Then 5573 rows should pass the filter", () => filterPasses(page, 5573));
    await session.step(34, "When user types \"RACE = Asian\" into table search", () => typeInto(page, "RACE = Asian", el("table search")));
    await session.step(35, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(36, "Then 72 rows should pass the filter", () => filterPasses(page, 72));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("Two conditions joined by \"and\" filter their intersection", {tag: ["@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(18, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await knownFailure(async () => {
      await session.step(42, "When user types \"AGE > 50 and SEX = M\" into table search", () => typeInto(page, "AGE > 50 and SEX = M", el("table search")));
      await session.step(43, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(44, "Then 856 rows should pass the filter", () => filterPasses(page, 856));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
  });
  test("Two conditions joined by \"or\" filter their union", {tag: ["@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(18, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await knownFailure(async () => {
      await session.step(50, "When user types \"AGE > 50 or SEX = M\" into table search", () => typeInto(page, "AGE > 50 or SEX = M", el("table search")));
      await session.step(51, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(52, "Then 3927 rows should pass the filter", () => filterPasses(page, 3927));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
  });
  test("Two conditions that cover every row keep every row", {tag: ["@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(18, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await knownFailure(async () => {
      await session.step(58, "When user types \"CONTROL = true\" into table search", () => typeInto(page, "CONTROL = true", el("table search")));
      await session.step(59, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(60, "Then 39 rows should pass the filter", () => filterPasses(page, 39));
      await session.step(61, "When user types \"SEX = M or SEX = F\" into table search", () => typeInto(page, "SEX = M or SEX = F", el("table search")));
      await session.step(62, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(63, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
  });
  test("A year alone compares a date column as the first day of that year does", {tag: ["@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(18, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await knownFailure(async () => {
      await session.step(70, "When user types \"STARTED > 1/1/1990\" into table search", () => typeInto(page, "STARTED > 1/1/1990", el("table search")));
      await session.step(71, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(72, "Then 5573 rows should pass the filter", () => filterPasses(page, 5573));
      await session.step(73, "When user clears table search", () => clearField(page, el("table search")));
      await session.step(74, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(75, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(76, "When user types \"STARTED > 1990\" into table search", () => typeInto(page, "STARTED > 1990", el("table search")));
      await session.step(77, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(78, "Then 5573 rows should pass the filter", () => filterPasses(page, 5573));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
    });
  });
  test("Text padded with spaces matches as the text itself", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(18, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await session.step(82, "When user types \"Asian\" into table search", () => typeInto(page, "Asian", el("table search")));
    await session.step(83, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(84, "Then 5339 rows should pass the filter", () => filterPasses(page, 5339));
    await session.step(85, "When user clears table search", () => clearField(page, el("table search")));
    await session.step(86, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(87, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(88, "When user types \" Asian\" into table search", () => typeInto(page, " Asian", el("table search")));
    await session.step(89, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(90, "Then 5339 rows should pass the filter", () => filterPasses(page, 5339));
    await session.step(91, "When user clears table search", () => clearField(page, el("table search")));
    await session.step(92, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(93, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(94, "When user types \"Asian \" into table search", () => typeInto(page, "Asian ", el("table search")));
    await session.step(95, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(96, "Then 5339 rows should pass the filter", () => filterPasses(page, 5339));
    await session.step(97, "When user clears table search", () => clearField(page, el("table search")));
    await session.step(98, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(99, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(100, "When user types \"  Asian  \" into table search", () => typeInto(page, "  Asian  ", el("table search")));
    await session.step(101, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(102, "Then 5339 rows should pass the filter", () => filterPasses(page, 5339));
    await session.step(103, "And no errors should have been logged", () => noErrors(page));
  });
  test("A quoted value matches as the unquoted one", {tag: ["@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(18, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await knownFailure(async () => {
      await session.step(108, "When user types \"RACE = Asian\" into table search", () => typeInto(page, "RACE = Asian", el("table search")));
      await session.step(109, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(110, "Then 72 rows should pass the filter", () => filterPasses(page, 72));
      await session.step(111, "When user clears table search", () => clearField(page, el("table search")));
      await session.step(112, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(113, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(114, "When user types 'RACE = \"Asian\"' into table search", () => typeInto(page, "RACE = \"Asian\"", el("table search")));
      await session.step(115, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
      await session.step(116, "Then 72 rows should pass the filter", () => filterPasses(page, 72));
      await session.step(117, "And no errors should have been logged", () => noErrors(page));
    });
  });
  test("A date condition over a datetime column with empty cells filters without an error", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(18, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await session.step(120, "Given user opens a table \"synth2019\" with:", () => openTableOf(page, "synth2019", [["EVENT_DATE","NUM"],["2018-12-15","1"],["","2"],["2019-03-01","3"],["2019-06-15","4"],["","5"],["2020-01-10","6"]]), [["EVENT_DATE","NUM"],["2018-12-15","1"],["","2"],["2019-03-01","3"],["2019-06-15","4"],["","5"],["2020-01-10","6"]]);
    await session.step(128, "Then \"EVENT_DATE\" column should have type \"datetime\"", () => columnType(page, "EVENT_DATE", "datetime"));
    await session.step(129, "And \"EVENT_DATE\" column should have missing values", () => columnIncomplete(page, "EVENT_DATE"));
    await session.step(130, "When user presses Control+f in grid overlay", () => pressKeyIn(page, "Control+f", el("grid overlay")));
    await session.step(131, "Then table search should be focused", () => shouldBe(page, el("table search"), "focused"));
    await session.step(132, "When user types \"NUM > 4\" into table search", () => typeInto(page, "NUM > 4", el("table search")));
    await session.step(133, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(134, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
    await session.step(135, "When user types \"EVENT_DATE > 1/1/2019\" into table search", () => typeInto(page, "EVENT_DATE > 1/1/2019", el("table search")));
    await session.step(136, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(137, "Then 3 rows should pass the filter", () => filterPasses(page, 3));
    await session.step(139, "When user types \"2019\" into table search", () => typeInto(page, "2019", el("table search")));
    await session.step(140, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(141, "Then 2 rows should pass the filter", () => filterPasses(page, 2));
    await session.step(142, "When user clears table search", () => clearField(page, el("table search")));
    await session.step(143, "And user presses Enter in table search", () => pressKeyIn(page, "Enter", el("table search")));
    await session.step(144, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(145, "And no errors should have been logged", () => noErrors(page));
  });
});
