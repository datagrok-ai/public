/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/data-text-to-columns.feature
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
import {check, clickOn, enterInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {cellIsCurrentObject, rowCount, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Split Column splits a text column into new columns by a delimiter", () => {
  const session = feature(test, "features/table-view/data-text-to-columns.feature", import.meta.url);
  test("Split by a delimiter adds one column per part, the shorter values leaving the rest empty", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens a table \"codes\" with:", () => openTableOf(page, "codes", [["code","n"],["a;b","1"],["c;d;e","2"],["f","3"]]), [["code","n"],["a;b","1"],["c;d;e","2"],["f","3"]]);
    await session.step(21, "And the \"code\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "code", 1));
    await session.step(24, "When user picks \"Data > Split Column...\" from the top menu", () => pickFromTopMenu(page, "Data > Split Column..."));
    await session.step(25, "Then Split Column dialog should be visible", () => shouldBe(page, el("Split Column dialog"), "visible"));
    await session.step(26, "And Column input in Split Column dialog should contain text \"code\"", () => shouldContainText(page, el("Column input in Split Column dialog"), "code"));
    await session.step(27, "When user enters \";\" into \"Split by\" input in Split Column dialog", () => enterInto(page, ";", el("\"Split by\" input in Split Column dialog")));
    await session.step(28, "And user clicks on OK button in Split Column dialog", () => clickOn(page, el("OK button in Split Column dialog")));
    await session.step(29, "Then Split Column dialog should be hidden", () => shouldBe(page, el("Split Column dialog"), "hidden"));
    await session.step(30, "And the table should have the columns \"code, n, col1, col2, col3\"", () => columnsExactly(page, "code, n, col1, col2, col3"));
    await session.step(31, "And the value of \"col1\" column in row 1 should be \"a\"", () => valueInRow(page, "col1", 1, "a"));
    await session.step(32, "And the value of \"col2\" column in row 1 should be \"b\"", () => valueInRow(page, "col2", 1, "b"));
    await session.step(33, "And the value of \"col3\" column in row 1 should be \"\"", () => valueInRow(page, "col3", 1, ""));
    await session.step(34, "And the value of \"col3\" column in row 2 should be \"e\"", () => valueInRow(page, "col3", 2, "e"));
    await session.step(35, "And the value of \"col1\" column in row 3 should be \"f\"", () => valueInRow(page, "col1", 3, "f"));
    await session.step(36, "And the value of \"col2\" column in row 3 should be \"\"", () => valueInRow(page, "col2", 3, ""));
    await session.step(37, "And the value of \"code\" column in row 2 should be \"c;d;e\"", () => valueInRow(page, "code", 2, "c;d;e"));
    await session.step(38, "And the table should have 3 rows", () => rowCount(page, 3));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("Prefix names the new columns", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens a table \"codes\" with:", () => openTableOf(page, "codes", [["code","n"],["a;b","1"],["c;d;e","2"],["f","3"]]), [["code","n"],["a;b","1"],["c;d;e","2"],["f","3"]]);
    await session.step(21, "And the \"code\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "code", 1));
    await session.step(42, "When user picks \"Data > Split Column...\" from the top menu", () => pickFromTopMenu(page, "Data > Split Column..."));
    await session.step(43, "And user enters \";\" into \"Split by\" input in Split Column dialog", () => enterInto(page, ";", el("\"Split by\" input in Split Column dialog")));
    await session.step(44, "And user enters \"part\" into \"Prefix\" input in Split Column dialog", () => enterInto(page, "part", el("\"Prefix\" input in Split Column dialog")));
    await session.step(45, "And user clicks on OK button in Split Column dialog", () => clickOn(page, el("OK button in Split Column dialog")));
    await session.step(46, "Then the table should have the columns \"code, n, part1, part2, part3\"", () => columnsExactly(page, "code, n, part1, part2, part3"));
    await session.step(47, "And the value of \"part2\" column in row 2 should be \"d\"", () => valueInRow(page, "part2", 2, "d"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
  test("By chars splits on each character of Split by", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens a table \"codes\" with:", () => openTableOf(page, "codes", [["code","n"],["a;b","1"],["c;d;e","2"],["f","3"]]), [["code","n"],["a;b","1"],["c;d;e","2"],["f","3"]]);
    await session.step(21, "And the \"code\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "code", 1));
    await session.step(51, "Given user sets \"code\" column in row 1 to \"a;b,c\"", () => setCell(page, "code", 1, "a;b,c"));
    await session.step(52, "When user picks \"Data > Split Column...\" from the top menu", () => pickFromTopMenu(page, "Data > Split Column..."));
    await session.step(53, "And user enters \";,\" into \"Split by\" input in Split Column dialog", () => enterInto(page, ";,", el("\"Split by\" input in Split Column dialog")));
    await session.step(54, "And user checks \"By chars\" input in Split Column dialog", () => check(page, el("\"By chars\" input in Split Column dialog")));
    await session.step(55, "And user clicks on OK button in Split Column dialog", () => clickOn(page, el("OK button in Split Column dialog")));
    await session.step(56, "Then the table should have the columns \"code, n, col1, col2, col3\"", () => columnsExactly(page, "code, n, col1, col2, col3"));
    await session.step(57, "And the value of \"col1\" column in row 1 should be \"a\"", () => valueInRow(page, "col1", 1, "a"));
    await session.step(58, "And the value of \"col2\" column in row 1 should be \"b\"", () => valueInRow(page, "col2", 1, "b"));
    await session.step(59, "And the value of \"col3\" column in row 1 should be \"c\"", () => valueInRow(page, "col3", 1, "c"));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
  });
  test("Consecutive as one treats repeated delimiters as one", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens a table \"codes\" with:", () => openTableOf(page, "codes", [["code","n"],["a;b","1"],["c;d;e","2"],["f","3"]]), [["code","n"],["a;b","1"],["c;d;e","2"],["f","3"]]);
    await session.step(21, "And the \"code\" cell of row 1 is the current object", () => cellIsCurrentObject(page, "code", 1));
    await session.step(63, "Given user sets \"code\" column in row 1 to \"a;;b\"", () => setCell(page, "code", 1, "a;;b"));
    await session.step(64, "When user picks \"Data > Split Column...\" from the top menu", () => pickFromTopMenu(page, "Data > Split Column..."));
    await session.step(65, "And user enters \";\" into \"Split by\" input in Split Column dialog", () => enterInto(page, ";", el("\"Split by\" input in Split Column dialog")));
    await session.step(66, "And user checks \"Consecutive as one\" input in Split Column dialog", () => check(page, el("\"Consecutive as one\" input in Split Column dialog")));
    await session.step(67, "And user clicks on OK button in Split Column dialog", () => clickOn(page, el("OK button in Split Column dialog")));
    await session.step(68, "Then the value of \"col1\" column in row 1 should be \"a\"", () => valueInRow(page, "col1", 1, "a"));
    await session.step(69, "And the value of \"col2\" column in row 1 should be \"b\"", () => valueInRow(page, "col2", 1, "b"));
    await session.step(70, "And the value of \"col3\" column in row 2 should be \"e\"", () => valueInRow(page, "col3", 2, "e"));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
  });
});
