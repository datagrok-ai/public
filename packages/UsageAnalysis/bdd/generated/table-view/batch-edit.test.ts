/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/batch-edit.feature
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
import {enterInto, selectIn, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {displayedInRowOfTable} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectInTableOneOf, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Batch Edit shows the values a set of rows shares and writes a typed value into every row of the set", () => {
  const session = feature(test, "features/table-view/batch-edit.feature", import.meta.url);
  test("The editor opens on the selected rows, shows their shared value, and a typed value goes into every selected row", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cid","Oslo","3"]]), [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cid","Oslo","3"]]);
    await session.step(26, "Given user selects rows where \"city\" is \"Oslo\"", () => selectWhereIs(page, "city", "Oslo"));
    await session.step(27, "When user picks \"Data > Batch Edit...\" from the top menu", () => pickFromTopMenu(page, "Data > Batch Edit..."));
    await session.step(28, "Then \"Rows\" input should have value \"Selected\"", () => shouldHaveValue(page, el("\"Rows\" input"), "Selected"));
    await session.step(29, "And value of \"city\" property should have text \"Oslo\"", () => shouldHaveText(page, el("value of \"city\" property"), "Oslo"));
    await session.step(30, "And value of \"name\" property should have text \"\"", () => shouldHaveText(page, el("value of \"name\" property"), ""));
    await session.step(31, "When user enters \"Paris\" into \"city\" property", () => enterInto(page, "Paris", el("\"city\" property")));
    await session.step(32, "Then the \"city\" cell of row 1 of table \"people\" should be displayed as \"Paris\"", () => displayedInRowOfTable(page, "city", 1, "people", "Paris"));
    await session.step(33, "And the \"city\" cell of row 3 of table \"people\" should be displayed as \"Paris\"", () => displayedInRowOfTable(page, "city", 3, "people", "Paris"));
    await session.step(34, "And the \"city\" cell of row 2 of table \"people\" should be displayed as \"Rome\"", () => displayedInRowOfTable(page, "city", 2, "people", "Rome"));
    await session.step(35, "And the \"name\" cell of row 1 of table \"people\" should be displayed as \"Ann\"", () => displayedInRowOfTable(page, "name", 1, "people", "Ann"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
  test("The editor follows the selection", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cid","Oslo","3"]]), [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cid","Oslo","3"]]);
    await session.step(39, "Given user selects rows where \"name\" is \"Ann\"", () => selectWhereIs(page, "name", "Ann"));
    await session.step(40, "When user picks \"Data > Batch Edit...\" from the top menu", () => pickFromTopMenu(page, "Data > Batch Edit..."));
    await session.step(41, "Then value of \"name\" property should have text \"Ann\"", () => shouldHaveText(page, el("value of \"name\" property"), "Ann"));
    await session.step(42, "When user selects rows of table \"people\" where \"name\" is one of \"Bob\"", () => selectInTableOneOf(page, "people", "name", "Bob"));
    await session.step(43, "Then value of \"name\" property should have text \"Bob\"", () => shouldHaveText(page, el("value of \"name\" property"), "Bob"));
    await session.step(44, "And value of \"city\" property should have text \"Rome\"", () => shouldHaveText(page, el("value of \"city\" property"), "Rome"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
  test("Rows set to All writes the value into every row", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cid","Oslo","3"]]), [["name","city","weight"],["Ann","Oslo","1.5"],["Bob","Rome","2.25"],["Cid","Oslo","3"]]);
    await session.step(48, "Given user selects rows where \"name\" is \"Ann\"", () => selectWhereIs(page, "name", "Ann"));
    await session.step(49, "When user picks \"Data > Batch Edit...\" from the top menu", () => pickFromTopMenu(page, "Data > Batch Edit..."));
    await session.step(50, "And user selects \"All\" in \"Rows\" input", () => selectIn(page, "All", el("\"Rows\" input")));
    await session.step(51, "Then value of \"city\" property should have text \"\"", () => shouldHaveText(page, el("value of \"city\" property"), ""));
    await session.step(52, "When user enters \"Paris\" into \"city\" property", () => enterInto(page, "Paris", el("\"city\" property")));
    await session.step(53, "Then the \"city\" cell of row 1 of table \"people\" should be displayed as \"Paris\"", () => displayedInRowOfTable(page, "city", 1, "people", "Paris"));
    await session.step(54, "And the \"city\" cell of row 2 of table \"people\" should be displayed as \"Paris\"", () => displayedInRowOfTable(page, "city", 2, "people", "Paris"));
    await session.step(55, "And the \"city\" cell of row 3 of table \"people\" should be displayed as \"Paris\"", () => displayedInRowOfTable(page, "city", 3, "people", "Paris"));
    await session.step(56, "And the \"name\" cell of row 2 of table \"people\" should be displayed as \"Bob\"", () => displayedInRowOfTable(page, "name", 2, "people", "Bob"));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
  });
});
