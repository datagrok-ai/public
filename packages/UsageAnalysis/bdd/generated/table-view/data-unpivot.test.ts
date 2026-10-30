/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/data-unpivot.feature
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
import {clickOn, enterInto, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount, tableColumns, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Unpivot turns the merged columns into category and value rows of a new table", () => {
  const session = feature(test, "features/table-view/data-unpivot.feature", import.meta.url);
  test("Copy and Merge give one row per student and subject in a new table view", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens a table \"scores\" with:", () => openTableOf(page, "scores", [["student","math","art"],["Ann","5","3"],["Bob","4","2"]]), [["student","math","art"],["Ann","5","3"],["Bob","4","2"]]);
    await session.step(23, "When user picks \"Data > Unpivot...\" from the top menu", () => pickFromTopMenu(page, "Data > Unpivot..."));
    await session.step(24, "Then Unpivot dialog should be visible", () => shouldBe(page, el("Unpivot dialog"), "visible"));
    await session.step(25, "And OK button in Unpivot dialog should be disabled", () => shouldBe(page, el("OK button in Unpivot dialog"), "disabled"));
    await session.step(26, "When user selects \"Copy\" in 3rd input in Unpivot dialog", () => selectIn(page, "Copy", el("3rd input in Unpivot dialog")));
    await session.step(27, "And user selects \"Merge\" in 4th input in Unpivot dialog", () => selectIn(page, "Merge", el("4th input in Unpivot dialog")));
    await session.step(28, "Then OK button in Unpivot dialog should be enabled", () => shouldBe(page, el("OK button in Unpivot dialog"), "enabled"));
    await session.step(29, "When user selects \"Merge\" in 5th input in Unpivot dialog", () => selectIn(page, "Merge", el("5th input in Unpivot dialog")));
    await session.step(30, "And user clicks on OK button in Unpivot dialog", () => clickOn(page, el("OK button in Unpivot dialog")));
    await session.step(31, "Then Unpivot dialog should be hidden", () => shouldBe(page, el("Unpivot dialog"), "hidden"));
    await session.step(32, "And the open table views should be exactly \"scores, result\"", () => openTableViewsExactly(page, "scores, result"));
    await session.step(33, "Given user switches to the \"result\" table view", () => switchTableView(page, "result"));
    await session.step(34, "Then the table should have 4 rows", () => rowCount(page, 4));
    await session.step(35, "And the table should have the columns \"student, Category, Value\"", () => columnsExactly(page, "student, Category, Value"));
    await session.step(36, "And the value of \"student\" column in row 1 should be \"Ann\"", () => valueInRow(page, "student", 1, "Ann"));
    await session.step(37, "And the value of \"Category\" column in row 1 should be \"math\"", () => valueInRow(page, "Category", 1, "math"));
    await session.step(38, "And the value of \"Value\" column in row 1 should be \"5\"", () => valueInRow(page, "Value", 1, "5"));
    await session.step(39, "And the value of \"Category\" column in row 2 should be \"art\"", () => valueInRow(page, "Category", 2, "art"));
    await session.step(40, "And the value of \"Value\" column in row 2 should be \"3\"", () => valueInRow(page, "Value", 2, "3"));
    await session.step(41, "And the value of \"student\" column in row 4 should be \"Bob\"", () => valueInRow(page, "student", 4, "Bob"));
    await session.step(42, "And the value of \"Value\" column in row 4 should be \"2\"", () => valueInRow(page, "Value", 4, "2"));
    await session.step(43, "And table \"scores\" should have 2 rows", () => tableRows(page, "scores", 2));
    await session.step(44, "And table \"scores\" should have columns \"student, math, art\"", () => tableColumns(page, "scores", "student, math, art"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
  test("The names typed in the dialog name the category and value columns, and a column left out is dropped", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens a table \"scores\" with:", () => openTableOf(page, "scores", [["student","math","art"],["Ann","5","3"],["Bob","4","2"]]), [["student","math","art"],["Ann","5","3"],["Bob","4","2"]]);
    await session.step(48, "When user picks \"Data > Unpivot...\" from the top menu", () => pickFromTopMenu(page, "Data > Unpivot..."));
    await session.step(49, "And user enters \"subject\" into \"Category column name\" input in Unpivot dialog", () => enterInto(page, "subject", el("\"Category column name\" input in Unpivot dialog")));
    await session.step(50, "And user enters \"mark\" into \"Value column name\" input in Unpivot dialog", () => enterInto(page, "mark", el("\"Value column name\" input in Unpivot dialog")));
    await session.step(51, "And user selects \"Merge\" in 4th input in Unpivot dialog", () => selectIn(page, "Merge", el("4th input in Unpivot dialog")));
    await session.step(52, "And user selects \"Merge\" in 5th input in Unpivot dialog", () => selectIn(page, "Merge", el("5th input in Unpivot dialog")));
    await session.step(53, "And user clicks on OK button in Unpivot dialog", () => clickOn(page, el("OK button in Unpivot dialog")));
    await session.step(54, "Then the open table views should be exactly \"scores, result\"", () => openTableViewsExactly(page, "scores, result"));
    await session.step(55, "Given user switches to the \"result\" table view", () => switchTableView(page, "result"));
    await session.step(56, "Then the table should have 4 rows", () => rowCount(page, 4));
    await session.step(57, "And the table should have the columns \"subject, mark\"", () => columnsExactly(page, "subject, mark"));
    await session.step(58, "And the value of \"subject\" column in row 2 should be \"art\"", () => valueInRow(page, "subject", 2, "art"));
    await session.step(59, "And the value of \"mark\" column in row 2 should be \"3\"", () => valueInRow(page, "mark", 2, "3"));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
});
