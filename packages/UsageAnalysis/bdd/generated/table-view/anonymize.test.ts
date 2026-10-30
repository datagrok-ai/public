/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/table-view/anonymize.feature
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
import {clickOn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueBetween, someValueDiffers, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {checkedInColumnList, toggleInColumnList, uncheckedInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Data > Anonymize replaces the values of the checked columns in place", () => {
  const session = feature(test, "features/table-view/anonymize.feature", import.meta.url);
  test("The dialog lists the columns unchecked, and OK replaces the categories of a checked text column only", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","w","w0"],["Ann","Oslo","10.5","10.5"],["Bob","Rome","10.5","10.5"],["Cid","Oslo","10.5","10.5"]]), [["name","city","w","w0"],["Ann","Oslo","10.5","10.5"],["Bob","Rome","10.5","10.5"],["Cid","Oslo","10.5","10.5"]]);
    await session.step(27, "When user picks \"Data > Anonymize...\" from the top menu", () => pickFromTopMenu(page, "Data > Anonymize..."));
    await session.step(28, "Then Anonymize Data dialog should be visible", () => shouldBe(page, el("Anonymize Data dialog"), "visible"));
    await session.step(29, "And \"Number randomization factor\" input in Anonymize Data dialog should have value \"0.1\"", () => shouldHaveValue(page, el("\"Number randomization factor\" input in Anonymize Data dialog"), "0.1"));
    await session.step(30, "And the \"city\" column should not be checked in the column list of Anonymize Data dialog", () => uncheckedInColumnList(page, "city", el("Anonymize Data dialog")));
    await session.step(31, "When user toggles the \"city\" column in the column list of Anonymize Data dialog", () => toggleInColumnList(page, "city", el("Anonymize Data dialog")));
    await session.step(32, "Then the \"city\" column should be checked in the column list of Anonymize Data dialog", () => checkedInColumnList(page, "city", el("Anonymize Data dialog")));
    await session.step(33, "When user clicks on OK button in Anonymize Data dialog", () => clickOn(page, el("OK button in Anonymize Data dialog")));
    await session.step(34, "Then Anonymize Data dialog should be hidden", () => shouldBe(page, el("Anonymize Data dialog"), "hidden"));
    await session.step(35, "And the value of \"city\" column in row 1 should be \"city 0\"", () => valueInRow(page, "city", 1, "city 0"));
    await session.step(36, "And the value of \"city\" column in row 2 should be \"city 1\"", () => valueInRow(page, "city", 2, "city 1"));
    await session.step(37, "And the value of \"city\" column in row 3 should be \"city 0\"", () => valueInRow(page, "city", 3, "city 0"));
    await session.step(38, "And the value of \"name\" column in row 1 should be \"Ann\"", () => valueInRow(page, "name", 1, "Ann"));
    await session.step(39, "And every value of \"w\" column should lie between 10.5 and 10.5", () => everyValueBetween(page, "w", 10.5, 10.5));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
  test("A checked number column is randomized within the factor's range, and an unchecked one is left alone", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And user opens a table \"people\" with:", () => openTableOf(page, "people", [["name","city","w","w0"],["Ann","Oslo","10.5","10.5"],["Bob","Rome","10.5","10.5"],["Cid","Oslo","10.5","10.5"]]), [["name","city","w","w0"],["Ann","Oslo","10.5","10.5"],["Bob","Rome","10.5","10.5"],["Cid","Oslo","10.5","10.5"]]);
    await session.step(43, "When user picks \"Data > Anonymize...\" from the top menu", () => pickFromTopMenu(page, "Data > Anonymize..."));
    await session.step(44, "And user toggles the \"w\" column in the column list of Anonymize Data dialog", () => toggleInColumnList(page, "w", el("Anonymize Data dialog")));
    await session.step(45, "And user clicks on OK button in Anonymize Data dialog", () => clickOn(page, el("OK button in Anonymize Data dialog")));
    await session.step(46, "Then Anonymize Data dialog should be hidden", () => shouldBe(page, el("Anonymize Data dialog"), "hidden"));
    await session.step(47, "And some value of \"w\" column should differ from \"w0\" column in the same row", () => someValueDiffers(page, "w", "w0"));
    await session.step(48, "And every value of \"w\" column should lie between 8.34 and 13.22", () => everyValueBetween(page, "w", 8.34, 13.22));
    await session.step(49, "And every value of \"w0\" column should lie between 10.5 and 10.5", () => everyValueBetween(page, "w0", 10.5, 10.5));
    await session.step(50, "And the value of \"city\" column in row 2 should be \"Rome\"", () => valueInRow(page, "city", 2, "Rome"));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
});
