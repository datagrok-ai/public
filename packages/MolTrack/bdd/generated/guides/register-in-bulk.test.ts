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
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveItemsInTotal} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: register compounds in bulk", () => {
  const session = feature(test, "features/guides/register-in-bulk.feature", import.meta.url);
  test("Map a file onto compounds, preview it, register it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the MolTrack EMS demo data is on the server", () => demoData(page));
    await session.step(15, "And the bulk registration example is open, and its compounds are not registered now and at feature end", () => bulkExample(page));
    await session.step(16, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"));
    enter(page, "MolTrack EMS");
    await session.step(18, "When user clicks on more actions in app ribbon", () => clickOn(page, el("more actions in app ribbon")), undefined, "Bulk registration is the Import of the Compounds view");
    await session.step(19, "And user clicks on \"Import…\" menu item", () => clickOn(page, el("\"Import…\" menu item")));
    await session.step(20, "Then \"Register compounds in bulk\" dialog should be visible", () => shouldBe(page, el("\"Register compounds in bulk\" dialog"), "visible"));
    await session.step(22, "And \"Data\" input in \"Register compounds in bulk\" dialog should contain text \"bulk-registration\"", () => shouldContainText(page, el("\"Data\" input in \"Register compounds in bulk\" dialog"), "bulk-registration"), undefined, "The table, and what a row with a problem does: here, it is skipped");
    await session.step(23, "When user selects \"Skip the row\" in \"On error\" input in \"Register compounds in bulk\" dialog", () => selectIn(page, "Skip the row", el("\"On error\" input in \"Register compounds in bulk\" dialog")));
    await session.step(24, "And user clicks on \"Next\" button in \"Register compounds in bulk\" dialog", () => clickOn(page, el("\"Next\" button in \"Register compounds in bulk\" dialog")));
    await session.step(26, "Then register mapping should contain text \"Common Name\"", () => shouldContainText(page, el("register mapping"), "Common Name"), undefined, "Each column of the file, and what it registers as");
    await session.step(27, "When user clicks on \"Next\" button in \"Register compounds in bulk\" dialog", () => clickOn(page, el("\"Next\" button in \"Register compounds in bulk\" dialog")));
    await session.step(29, "Then register preview should contain text \"6 rows: 3 new compounds, 1 repeating an earlier row, 1 registered already, 1 with a problem\"", () => shouldContainText(page, el("register preview"), "6 rows: 3 new compounds, 1 repeating an earlier row, 1 registered already, 1 with a problem"), undefined, "The preview: each row by its name, and what it becomes");
    await session.step(30, "And register preview should contain text \"Sodium hentriacontanoate\"", () => shouldContainText(page, el("register preview"), "Sodium hentriacontanoate"));
    await session.step(31, "And register preview should contain text \"Same compound as row 1\"", () => shouldContainText(page, el("register preview"), "Same compound as row 1"));
    await session.step(32, "And register preview should contain text \"This structure is already registered as DG-\"", () => shouldContainText(page, el("register preview"), "This structure is already registered as DG-"));
    await session.step(33, "And register preview should contain text \"Cannot read the structure\"", () => shouldContainText(page, el("register preview"), "Cannot read the structure"));
    await session.step(34, "And register preview should contain text \"The rows with a problem are skipped.\"", () => shouldContainText(page, el("register preview"), "The rows with a problem are skipped."));
    await session.step(35, "When user clicks on \"Register\" button in \"Register compounds in bulk\" dialog", () => clickOn(page, el("\"Register\" button in \"Register compounds in bulk\" dialog")));
    await session.step(37, "Then register report should contain text \"Registered 3 of 6 rows\"", () => shouldContainText(page, el("register report"), "Registered 3 of 6 rows"), undefined, "The report, row by row, with each row's Grok ID");
    await session.step(38, "And register report should contain text \"Same compound as row 1\"", () => shouldContainText(page, el("register report"), "Same compound as row 1"));
    await session.step(39, "And register report should contain text \"DG-\"", () => shouldContainText(page, el("register report"), "DG-"));
    await session.step(41, "When user clicks on \"Show the 3 compounds registered\" link in register report", () => clickOn(page, el("\"Show the 3 compounds registered\" link in register report")), undefined, "Open the compounds just registered");
    await session.step(42, "Then compound list should have 3 items in total", () => shouldHaveItemsInTotal(page, el("compound list"), 3));
    await session.step(43, "And compound list should contain text \"Hentriacontanoic acid\"", () => shouldContainText(page, el("compound list"), "Hentriacontanoic acid"));
    await session.step(44, "And the bulk registration example should have 3 new compounds on the server", () => bulkRegistered(page, 3));
  });
});
