/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/register-in-bulk.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bulkExample, bulkRegistered, demoData, openEmsView} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldNotContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTypeOfTable} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: register compounds in bulk", () => {
  const session = feature(test, "features/guides/register-in-bulk.feature", import.meta.url);
  test("Pick a table, see what is new, register it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page), undefined, "Sign in");
    await session.step(13, "And simple mode is off", () => simpleModeOff(page), undefined, "Simple mode is off");
    await session.step(15, "And the MolTrack EMS demo data is on the server", () => demoData(page), undefined, "The demo data is loaded");
    await session.step(17, "And the bulk registration example is open, and its compounds are not registered now and at feature end", () => bulkExample(page), undefined, "Open the example table: three new structures, one already registered, one broken");
    await session.step(19, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"), undefined, "Open the Compounds view of MolTrack");
    enter(page, "MolTrack EMS");
    await session.step(21, "When user clicks on more actions in app ribbon", () => clickOn(page, el("more actions in app ribbon")), undefined, "Open the more actions menu");
    await session.step(23, "And user clicks on \"Import…\" menu item", () => clickOn(page, el("\"Import…\" menu item")), undefined, "Choose Import");
    await session.step(25, "Then \"Import compounds\" dialog should be visible", () => shouldBe(page, el("\"Import compounds\" dialog"), "visible"), undefined, "The import wizard opens");
    await session.step(27, "When user selects \"bulk-registration\" in \"Data\" input in \"Import compounds\" dialog", () => selectIn(page, "bulk-registration", el("\"Data\" input in \"Import compounds\" dialog")), undefined, "Pick the table to register");
    await session.step(29, "And user unchecks \"All or nothing\" input in \"Import compounds\" dialog", () => uncheck(page, el("\"All or nothing\" input in \"Import compounds\" dialog")), undefined, "Turn off all or nothing, so one broken row does not stop the rest");
    await session.step(31, "And user clicks on \"Next\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Next\" button in \"Import compounds\" dialog")), undefined, "Next");
    await session.step(33, "Then import mapping should contain text \"smiles\"", () => shouldContainText(page, el("import mapping"), "smiles"), undefined, "The smiles column fills the structure");
    await session.step(35, "And import mapping should not contain text \"Registration Hash\"", () => shouldNotContainText(page, el("import mapping"), "Registration Hash"), undefined, "No hash or matching rule to map: they are worked out");
    await session.step(37, "When user clicks on \"Next\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Next\" button in \"Import compounds\" dialog")), undefined, "Next");
    await session.step(39, "Then import preview should contain text \"3 will be added, 0 updated, 1 skipped, 1 row has errors\"", () => shouldContainText(page, el("import preview"), "3 will be added, 0 updated, 1 skipped, 1 row has errors"), undefined, "The preview: three will be added, one skipped, one with an error");
    await session.step(41, "And import preview should contain text \"Already registered as DG-\"", () => shouldContainText(page, el("import preview"), "Already registered as DG-"), undefined, "Row four is already registered, with its Grok ID");
    await session.step(43, "And import preview should contain text \"Cannot read the structure\"", () => shouldContainText(page, el("import preview"), "Cannot read the structure"), undefined, "Row five cannot be read");
    await session.step(45, "When user clicks on \"Import 3 compounds\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Import 3 compounds\" button in \"Import compounds\" dialog")), undefined, "Import the three new compounds");
    await session.step(47, "Then import report should contain text \"3 inserted, 0 updated, 1 skipped, 1 failed\"", () => shouldContainText(page, el("import report"), "3 inserted, 0 updated, 1 skipped, 1 failed"), undefined, "The report: three inserted, one skipped, one failed");
    await session.step(49, "And import report should contain text \"Grok ID\"", () => shouldContainText(page, el("import report"), "Grok ID"), undefined, "Every row has its Grok ID or the reason it has none");
    await session.step(51, "When user clicks on \"Close\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Close\" button in \"Import compounds\" dialog")), undefined, "Close the wizard");
    await session.step(53, "Then \"Grok ID\" column of table \"bulk-registration\" should have type \"string\"", () => columnTypeOfTable(page, "Grok ID", "bulk-registration", "string"), undefined, "The Grok IDs are written back to the table");
    await session.step(55, "And the bulk registration example should have 3 new compounds on the server", () => bulkRegistered(page, 3), undefined, "The three new compounds are on the server");
  });
});
