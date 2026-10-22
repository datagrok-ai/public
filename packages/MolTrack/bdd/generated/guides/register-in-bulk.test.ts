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
import {clickOn, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveItemsInTotal} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: register compounds in bulk", () => {
  const session = feature(test, "features/guides/register-in-bulk.feature", import.meta.url);
  test("Map a file onto compounds, preview it, register it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And the MolTrack EMS demo data is on the server", () => demoData(page));
    await session.step(16, "And the bulk registration example is open, and its compounds are not registered now and at feature end", () => bulkExample(page));
    await session.step(17, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"));
    enter(page, "MolTrack EMS");
    await session.step(19, "When user clicks on more actions in app ribbon", () => clickOn(page, el("more actions in app ribbon")), undefined, "Bulk registration is the Import of the Compounds view");
    await session.step(20, "And user clicks on \"Import…\" menu item", () => clickOn(page, el("\"Import…\" menu item")));
    await session.step(21, "Then \"Import compounds\" dialog should be visible", () => shouldBe(page, el("\"Import compounds\" dialog"), "visible"));
    await session.step(23, "And \"Data\" input in \"Import compounds\" dialog should contain text \"bulk-registration\"", () => shouldContainText(page, el("\"Data\" input in \"Import compounds\" dialog"), "bulk-registration"), undefined, "The table, and what a row with an error does: with \"All or nothing\" off, it fails alone");
    await session.step(24, "When user unchecks \"All or nothing\" input in \"Import compounds\" dialog", () => uncheck(page, el("\"All or nothing\" input in \"Import compounds\" dialog")));
    await session.step(25, "And user clicks on \"Next\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Next\" button in \"Import compounds\" dialog")));
    await session.step(27, "Then import mapping should contain text \"Common Name\"", () => shouldContainText(page, el("import mapping"), "Common Name"), undefined, "Each column of the file, and the compound column it fills");
    await session.step(28, "When user clicks on \"Next\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Next\" button in \"Import compounds\" dialog")));
    await session.step(30, "Then import preview should contain text \"3 will be added, 0 updated, 1 skipped, 1 row has errors.\"", () => shouldContainText(page, el("import preview"), "3 will be added, 0 updated, 1 skipped, 1 row has errors."), undefined, "The preview: what each row becomes");
    await session.step(31, "And import preview should contain text \"This structure is already registered as DG-\"", () => shouldContainText(page, el("import preview"), "This structure is already registered as DG-"));
    await session.step(32, "And import preview should contain text \"Cannot read the structure\"", () => shouldContainText(page, el("import preview"), "Cannot read the structure"));
    await session.step(33, "When user clicks on \"Next\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Next\" button in \"Import compounds\" dialog")));
    await session.step(35, "Then import report should contain text \"3 inserted, 0 updated, 1 skipped, 1 failed.\"", () => shouldContainText(page, el("import report"), "3 inserted, 0 updated, 1 skipped, 1 failed."), undefined, "The report, with each new compound's Grok ID");
    await session.step(36, "And import report should contain text \"DG-\"", () => shouldContainText(page, el("import report"), "DG-"));
    await session.step(38, "When user clicks on \"Show the 3 imported rows\" link in import report", () => clickOn(page, el("\"Show the 3 imported rows\" link in import report")), undefined, "Show the compounds just registered");
    await session.step(39, "Then compound list should have 3 items in total", () => shouldHaveItemsInTotal(page, el("compound list"), 3));
    await session.step(40, "And compound list should contain text \"Hentriacontanoic acid\"", () => shouldContainText(page, el("compound list"), "Hentriacontanoic acid"));
    await session.step(41, "And the bulk registration example should have 3 new compounds on the server", () => bulkRegistered(page, 3));
  });
});
