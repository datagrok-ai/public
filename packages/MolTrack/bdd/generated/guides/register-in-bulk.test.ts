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

test.describe("MolTrack EMS: register compounds and batches in bulk", () => {
  const session = feature(test, "features/guides/register-in-bulk.feature", import.meta.url);
  test("Map a file onto compounds and batches, preview it, register it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the MolTrack EMS demo data is on the server", () => demoData(page));
    await session.step(14, "And the bulk registration example is open, and its compounds are not registered now and at feature end", () => bulkExample(page));
    await session.step(15, "And user opens the \"Batches\" view of MolTrack EMS", () => openEmsView(page, "Batches"));
    enter(page, "MolTrack EMS");
    await session.step(17, "When user clicks on more actions in app ribbon", () => clickOn(page, el("more actions in app ribbon")), undefined, "Bulk registration is the Import of the Batches view");
    await session.step(18, "And user clicks on \"Import…\" menu item", () => clickOn(page, el("\"Import…\" menu item")));
    await session.step(19, "Then \"Register in bulk\" dialog should be visible", () => shouldBe(page, el("\"Register in bulk\" dialog"), "visible"));
    await session.step(21, "And \"Data\" input in \"Register in bulk\" dialog should contain text \"bulk-registration\"", () => shouldContainText(page, el("\"Data\" input in \"Register in bulk\" dialog"), "bulk-registration"), undefined, "The table, what to register, and what a row with a problem does: here, it is skipped");
    await session.step(22, "When user selects \"Skip the row\" in \"On error\" input in \"Register in bulk\" dialog", () => selectIn(page, "Skip the row", el("\"On error\" input in \"Register in bulk\" dialog")));
    await session.step(23, "And user clicks on \"Next\" button in \"Register in bulk\" dialog", () => clickOn(page, el("\"Next\" button in \"Register in bulk\" dialog")));
    await session.step(25, "Then register mapping should contain text \"Compound · Common Name\"", () => shouldContainText(page, el("register mapping"), "Compound · Common Name"), undefined, "Each column of the file, and what it registers as");
    await session.step(26, "And register mapping should contain text \"Batch · EPA Batch ID\"", () => shouldContainText(page, el("register mapping"), "Batch · EPA Batch ID"));
    await session.step(27, "When user clicks on \"Next\" button in \"Register in bulk\" dialog", () => clickOn(page, el("\"Next\" button in \"Register in bulk\" dialog")));
    await session.step(29, "Then register preview should contain text \"6 rows: 3 new compounds, 5 batches (1 of a registered compound), 1 with a problem\"", () => shouldContainText(page, el("register preview"), "6 rows: 3 new compounds, 5 batches (1 of a registered compound), 1 with a problem"), undefined, "The preview: each row by its name and batch ID, and what it becomes");
    await session.step(30, "And register preview should contain text \"Sodium hentriacontanoate\"", () => shouldContainText(page, el("register preview"), "Sodium hentriacontanoate"));
    await session.step(31, "And register preview should contain text \"Same compound as row 1\"", () => shouldContainText(page, el("register preview"), "Same compound as row 1"));
    await session.step(32, "And register preview should contain text \"Registered as DG-\"", () => shouldContainText(page, el("register preview"), "Registered as DG-"));
    await session.step(33, "And register preview should contain text \"Cannot read the structure\"", () => shouldContainText(page, el("register preview"), "Cannot read the structure"));
    await session.step(34, "And register preview should contain text \"The rows with a problem are skipped.\"", () => shouldContainText(page, el("register preview"), "The rows with a problem are skipped."));
    await session.step(35, "When user clicks on \"Register\" button in \"Register in bulk\" dialog", () => clickOn(page, el("\"Register\" button in \"Register in bulk\" dialog")));
    await session.step(37, "Then register report should contain text \"Registered 5 of 6 rows\"", () => shouldContainText(page, el("register report"), "Registered 5 of 6 rows"), undefined, "The report, row by row, with the Grok IDs given");
    await session.step(38, "And register report should contain text \"BDD-004\"", () => shouldContainText(page, el("register report"), "BDD-004"));
    await session.step(39, "And register report should contain text \"DGB-\"", () => shouldContainText(page, el("register report"), "DGB-"));
    await session.step(41, "When user clicks on \"Show the 5 batches registered\" link in register report", () => clickOn(page, el("\"Show the 5 batches registered\" link in register report")), undefined, "Open the batches just registered");
    await session.step(42, "Then batch list should have 5 items in total", () => shouldHaveItemsInTotal(page, el("batch list"), 5));
    await session.step(43, "And batch list should contain text \"BDD-001\"", () => shouldContainText(page, el("batch list"), "BDD-001"));
    await session.step(44, "And the bulk registration example should have 3 new compounds and 5 batches on the server", () => bulkRegistered(page, 3, 5));
  });
});
