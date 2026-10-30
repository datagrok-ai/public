/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/compare-tables.feature
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
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu, topMenuListsDisabled, topMenuListsEnabled} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {currentRowOfTableValue, rowCount, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {areaColor, areaNotColor, areasDiffer, clickArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Compare Tables opens the rows two tables differ in as a new color-coded \"diff\" table view", () => {
  const session = feature(test, "features/table-view/compare-tables.feature", import.meta.url);
  test("Compare Tables needs a second open table", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]);
    await session.step(28, "Then the top menu should list \"Data > Compare Tables...\" as disabled", () => topMenuListsDisabled(page, "Data > Compare Tables..."));
    await session.step(29, "Given user opens a table \"people2\" with:", () => openTableOf(page, "people2", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]);
    await session.step(34, "Then the top menu should list \"Data > Compare Tables...\" as enabled", () => topMenuListsEnabled(page, "Data > Compare Tables..."));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
  });
  test("The dialog opens keyed by name, and OK opens the rows that differ as a new table view", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]);
    await session.step(38, "Given user opens a table \"people2\" with:", () => openTableOf(page, "people2", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]);
    await session.step(43, "When user picks \"Data > Compare Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Tables..."));
    await session.step(44, "Then Compare Tables dialog should be visible", () => shouldBe(page, el("Compare Tables dialog"), "visible"));
    await session.step(45, "And Key Columns input in Compare Tables dialog should contain text \"name\"", () => shouldContainText(page, el("Key Columns input in Compare Tables dialog"), "name"));
    await session.step(46, "And OK button in Compare Tables dialog should be enabled", () => shouldBe(page, el("OK button in Compare Tables dialog"), "enabled"));
    await session.step(47, "When user clicks on OK button in Compare Tables dialog", () => clickOn(page, el("OK button in Compare Tables dialog")));
    await session.step(48, "Then Compare Tables dialog should be hidden", () => shouldBe(page, el("Compare Tables dialog"), "hidden"));
    await session.step(49, "And the open table views should be exactly \"people, people2, diff\"", () => openTableViewsExactly(page, "people, people2, diff"));
    await session.step(50, "Given user switches to the \"diff\" table view", () => switchTableView(page, "diff"));
    await session.step(51, "Then the table should have 3 rows", () => rowCount(page, 3));
    await session.step(52, "And the table should have the columns \"name, city, weight, ~t1 row, ~t2 row\"", () => columnsExactly(page, "name, city, weight, ~t1 row, ~t2 row"));
    await session.step(53, "And the value of \"name\" column in row 1 should be \"Bob\"", () => valueInRow(page, "name", 1, "Bob"));
    await session.step(54, "And the value of \"city\" column in row 1 should be \"Paris\"", () => valueInRow(page, "city", 1, "Paris"));
    await session.step(55, "And the value of \"name\" column in row 2 should be \"Cid\"", () => valueInRow(page, "name", 2, "Cid"));
    await session.step(56, "And the value of \"name\" column in row 3 should be \"Dan\"", () => valueInRow(page, "name", 3, "Dan"));
    await session.step(57, "And the value of \"weight\" column in row 3 should be \"4\"", () => valueInRow(page, "weight", 3, "4"));
    await session.step(58, "And table \"people\" should have 3 rows", () => tableRows(page, "people", 3));
    await session.step(59, "And table \"people2\" should have 3 rows", () => tableRows(page, "people2", 3));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
  test("The diff grid paints a changed value orange and the one-table rows in two other colors, and a click on a row makes it current in both tables", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]);
    await session.step(63, "Given user opens a table \"people2\" with:", () => openTableOf(page, "people2", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]);
    await session.step(68, "When user picks \"Data > Compare Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Tables..."));
    await session.step(69, "And user clicks on OK button in Compare Tables dialog", () => clickOn(page, el("OK button in Compare Tables dialog")));
    await session.step(70, "Then the open table views should be exactly \"people, people2, diff\"", () => openTableViewsExactly(page, "people, people2, diff"));
    await session.step(71, "Given user switches to the \"diff\" table view", () => switchTableView(page, "diff"));
    await session.step(72, "Then the \"cell 1 of city\" area of grid should contain the color \"#FFA500\"", () => areaColor(page, "cell 1 of city", el("grid"), "#FFA500"));
    await session.step(73, "And the \"cell 1 of name\" area of grid should not contain the color \"#FFA500\"", () => areaNotColor(page, "cell 1 of name", el("grid"), "#FFA500"));
    await session.step(74, "And the \"cell 2 of name\" and \"cell 3 of name\" areas of grid should be painted in different colors", () => areasDiffer(page, "cell 2 of name", "cell 3 of name", el("grid")));
    await session.step(75, "When user clicks on the \"cell 1 of name\" area of grid", () => clickArea(page, "cell 1 of name", el("grid")));
    await session.step(76, "Then \"name\" of the current row of table \"people\" should be \"Bob\"", () => currentRowOfTableValue(page, "name", "people", "Bob"));
    await session.step(77, "And \"name\" of the current row of table \"people2\" should be \"Bob\"", () => currentRowOfTableValue(page, "name", "people2", "Bob"));
    await session.step(78, "And no errors should have been logged", () => noErrors(page));
  });
  test("A key whose values repeat disables OK and says why", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]);
    await session.step(81, "Given user opens a table \"people2\" with:", () => openTableOf(page, "people2", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]);
    await session.step(86, "When user picks \"Data > Compare Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Tables..."));
    await session.step(87, "And user selects \"city\" in Key Columns input in Compare Tables dialog", () => selectIn(page, "city", el("Key Columns input in Compare Tables dialog")));
    await session.step(88, "Then Compare Tables dialog should contain text \"city is not unique\"", () => shouldContainText(page, el("Compare Tables dialog"), "city is not unique"));
    await session.step(89, "And OK button in Compare Tables dialog should be disabled", () => shouldBe(page, el("OK button in Compare Tables dialog"), "disabled"));
    await session.step(90, "When user clicks on CANCEL button in Compare Tables dialog", () => clickOn(page, el("CANCEL button in Compare Tables dialog")));
    await session.step(91, "Then the open table views should be exactly \"people, people2\"", () => openTableViewsExactly(page, "people, people2"));
    await session.step(92, "And no errors should have been logged", () => noErrors(page));
  });
  test("With no value columns only the rows in one table go into the diff", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Rome","10.5"],["Cid","Oslo","10.5"]]);
    await session.step(95, "Given user opens a table \"people2\" with:", () => openTableOf(page, "people2", [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]), [["name","city","weight"],["Ann","Oslo","10.5"],["Bob","Paris","10.5"],["Dan","Oslo","4"]]);
    await session.step(100, "When user picks \"Data > Compare Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Compare Tables..."));
    await session.step(101, "And user clicks on \"Clear\" icon in Compare Tables dialog", () => clickOn(page, el("\"Clear\" icon in Compare Tables dialog")));
    await session.step(102, "And user clicks on OK button in Compare Tables dialog", () => clickOn(page, el("OK button in Compare Tables dialog")));
    await session.step(103, "Then the open table views should be exactly \"people, people2, diff\"", () => openTableViewsExactly(page, "people, people2, diff"));
    await session.step(104, "Given user switches to the \"diff\" table view", () => switchTableView(page, "diff"));
    await session.step(105, "Then the table should have 2 rows", () => rowCount(page, 2));
    await session.step(106, "And the value of \"name\" column in row 1 should be \"Cid\"", () => valueInRow(page, "name", 1, "Cid"));
    await session.step(107, "And the value of \"name\" column in row 2 should be \"Dan\"", () => valueInRow(page, "name", 2, "Dan"));
    await session.step(108, "And no errors should have been logged", () => noErrors(page));
  });
});
