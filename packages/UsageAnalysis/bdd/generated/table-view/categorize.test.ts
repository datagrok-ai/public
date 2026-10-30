/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/categorize.feature
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
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, columnType, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {colorCodedAs, colorCodedCategorically, noColorCoding} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Categorize bins a number column into categories and saves them as a column, a matrix or a color coding", () => {
  const session = feature(test, "features/table-view/categorize.feature", import.meta.url);
  test("Bin fills the patterns, and SAVE COLUMN adds a color-coded category column", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","v"],["Ann","Oslo","1"],["Bob","Rome","2"],["Cid","Oslo","3"],["Dan","Rome","4"]]), [["name","city","v"],["Ann","Oslo","1"],["Bob","Rome","2"],["Cid","Oslo","3"],["Dan","Rome","4"]]);
    await session.step(28, "When user picks \"Data > Categorize...\" from the top menu", () => pickFromTopMenu(page, "Data > Categorize..."));
    await session.step(29, "Then Categorize Data dialog should be visible", () => shouldBe(page, el("Categorize Data dialog"), "visible"));
    await session.step(30, "When user selects \"v\" in \"Categorize\" input in Categorize Data dialog", () => selectIn(page, "v", el("\"Categorize\" input in Categorize Data dialog")));
    await session.step(31, "Then Categorize Data dialog should contain text \"Matching: 0\"", () => shouldContainText(page, el("Categorize Data dialog"), "Matching: 0"));
    await session.step(32, "When user clicks on \"Bin...\" link in Categorize Data dialog", () => clickOn(page, el("\"Bin...\" link in Categorize Data dialog")));
    await session.step(33, "Then OK button in Edit dialog should be disabled", () => shouldBe(page, el("OK button in Edit dialog"), "disabled"));
    await session.step(34, "When user enters \"2\" into \"Number of bins\" input in Edit dialog", () => enterInto(page, "2", el("\"Number of bins\" input in Edit dialog")));
    await session.step(35, "Then OK button in Edit dialog should be enabled", () => shouldBe(page, el("OK button in Edit dialog"), "enabled"));
    await session.step(36, "When user clicks on OK button in Edit dialog", () => clickOn(page, el("OK button in Edit dialog")));
    await session.step(37, "Then Edit dialog should be hidden", () => shouldBe(page, el("Edit dialog"), "hidden"));
    await session.step(38, "And Categorize Data dialog should contain text \"Matching: 4\"", () => shouldContainText(page, el("Categorize Data dialog"), "Matching: 4"));
    await session.step(39, "And Categorize Data dialog should contain text \"Not matching: 0\"", () => shouldContainText(page, el("Categorize Data dialog"), "Not matching: 0"));
    await session.step(40, "When user clicks on \"SAVE COLUMN\" button in Categorize Data dialog", () => clickOn(page, el("\"SAVE COLUMN\" button in Categorize Data dialog")));
    await session.step(41, "Then the value of \"category\" column in row 1 should be \"1 - 2.50\"", () => valueInRow(page, "category", 1, "1 - 2.50"));
    await session.step(42, "And the value of \"category\" column in row 2 should be \"1 - 2.50\"", () => valueInRow(page, "category", 2, "1 - 2.50"));
    await session.step(43, "And the value of \"category\" column in row 3 should be \"2.50 - 4\"", () => valueInRow(page, "category", 3, "2.50 - 4"));
    await session.step(44, "And the value of \"category\" column in row 4 should be \"2.50 - 4\"", () => valueInRow(page, "category", 4, "2.50 - 4"));
    await session.step(45, "And \"category\" column should be color-coded categorically", () => colorCodedCategorically(page, "category"));
    await session.step(46, "And \"v\" column should have no color coding", () => noColorCoding(page, "v"));
    await session.step(47, "And Categorize Data dialog should be visible", () => shouldBe(page, el("Categorize Data dialog"), "visible"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
  test("With Color Coding off, SAVE COLUMN adds the category column without colors", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","v"],["Ann","Oslo","1"],["Bob","Rome","2"],["Cid","Oslo","3"],["Dan","Rome","4"]]), [["name","city","v"],["Ann","Oslo","1"],["Bob","Rome","2"],["Cid","Oslo","3"],["Dan","Rome","4"]]);
    await session.step(51, "When user picks \"Data > Categorize...\" from the top menu", () => pickFromTopMenu(page, "Data > Categorize..."));
    await session.step(52, "And user selects \"v\" in \"Categorize\" input in Categorize Data dialog", () => selectIn(page, "v", el("\"Categorize\" input in Categorize Data dialog")));
    await session.step(53, "And user clicks on \"Bin...\" link in Categorize Data dialog", () => clickOn(page, el("\"Bin...\" link in Categorize Data dialog")));
    await session.step(54, "And user enters \"2\" into \"Number of bins\" input in Edit dialog", () => enterInto(page, "2", el("\"Number of bins\" input in Edit dialog")));
    await session.step(55, "And user clicks on OK button in Edit dialog", () => clickOn(page, el("OK button in Edit dialog")));
    await session.step(56, "And user unchecks \"Color Coding\" input in Categorize Data dialog", () => uncheck(page, el("\"Color Coding\" input in Categorize Data dialog")));
    await session.step(57, "And user clicks on \"SAVE COLUMN\" button in Categorize Data dialog", () => clickOn(page, el("\"SAVE COLUMN\" button in Categorize Data dialog")));
    await session.step(58, "Then the value of \"category\" column in row 4 should be \"2.50 - 4\"", () => valueInRow(page, "category", 4, "2.50 - 4"));
    await session.step(59, "And \"category\" column should have no color coding", () => noColorCoding(page, "category"));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
  test("SAVE MATRIX adds a boolean column per category", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","v"],["Ann","Oslo","1"],["Bob","Rome","2"],["Cid","Oslo","3"],["Dan","Rome","4"]]), [["name","city","v"],["Ann","Oslo","1"],["Bob","Rome","2"],["Cid","Oslo","3"],["Dan","Rome","4"]]);
    await session.step(63, "When user picks \"Data > Categorize...\" from the top menu", () => pickFromTopMenu(page, "Data > Categorize..."));
    await session.step(64, "And user selects \"v\" in \"Categorize\" input in Categorize Data dialog", () => selectIn(page, "v", el("\"Categorize\" input in Categorize Data dialog")));
    await session.step(65, "And user clicks on \"Bin...\" link in Categorize Data dialog", () => clickOn(page, el("\"Bin...\" link in Categorize Data dialog")));
    await session.step(66, "And user enters \"2\" into \"Number of bins\" input in Edit dialog", () => enterInto(page, "2", el("\"Number of bins\" input in Edit dialog")));
    await session.step(67, "And user clicks on OK button in Edit dialog", () => clickOn(page, el("OK button in Edit dialog")));
    await session.step(68, "And user clicks on \"SAVE MATRIX\" button in Categorize Data dialog", () => clickOn(page, el("\"SAVE MATRIX\" button in Categorize Data dialog")));
    await session.step(69, "Then \"1 - 2.50\" column should have type \"bool\"", () => columnType(page, "1 - 2.50", "bool"));
    await session.step(70, "And the value of \"1 - 2.50\" column in row 1 should be \"true\"", () => valueInRow(page, "1 - 2.50", 1, "true"));
    await session.step(71, "And the value of \"1 - 2.50\" column in row 3 should be \"false\"", () => valueInRow(page, "1 - 2.50", 3, "false"));
    await session.step(72, "And the value of \"2.50 - 4\" column in row 1 should be \"false\"", () => valueInRow(page, "2.50 - 4", 1, "false"));
    await session.step(73, "And the value of \"2.50 - 4\" column in row 4 should be \"true\"", () => valueInRow(page, "2.50 - 4", 4, "true"));
    await session.step(74, "And the table should have 5 columns", () => columnCount(page, 5));
    await session.step(75, "And no errors should have been logged", () => noErrors(page));
  });
  test("APPLY COLOR CODING colors the source column by the categories", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","v"],["Ann","Oslo","1"],["Bob","Rome","2"],["Cid","Oslo","3"],["Dan","Rome","4"]]), [["name","city","v"],["Ann","Oslo","1"],["Bob","Rome","2"],["Cid","Oslo","3"],["Dan","Rome","4"]]);
    await session.step(78, "When user picks \"Data > Categorize...\" from the top menu", () => pickFromTopMenu(page, "Data > Categorize..."));
    await session.step(79, "And user selects \"v\" in \"Categorize\" input in Categorize Data dialog", () => selectIn(page, "v", el("\"Categorize\" input in Categorize Data dialog")));
    await session.step(80, "And user clicks on \"Bin...\" link in Categorize Data dialog", () => clickOn(page, el("\"Bin...\" link in Categorize Data dialog")));
    await session.step(81, "And user enters \"2\" into \"Number of bins\" input in Edit dialog", () => enterInto(page, "2", el("\"Number of bins\" input in Edit dialog")));
    await session.step(82, "And user clicks on OK button in Edit dialog", () => clickOn(page, el("OK button in Edit dialog")));
    await session.step(83, "And user clicks on \"APPLY COLOR CODING\" button in Categorize Data dialog", () => clickOn(page, el("\"APPLY COLOR CODING\" button in Categorize Data dialog")));
    await session.step(84, "Then \"v\" column should be color-coded conditionally", () => colorCodedAs(page, "v", "conditionally"));
    await session.step(85, "And the table should have 3 columns", () => columnCount(page, 3));
    await session.step(86, "And no errors should have been logged", () => noErrors(page));
  });
});
