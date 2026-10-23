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
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTypeOfTable} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: register compounds in bulk", () => {
  const session = feature(test, "features/guides/register-in-bulk.feature", import.meta.url);
  test("Pick a table, see what is new, register it", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the MolTrack EMS demo data is on the server", () => demoData(page));
    await session.step(13, "And the bulk registration example is open, and its compounds are not registered now and at feature end", () => bulkExample(page));
    await session.step(14, "And user opens the \"Compounds\" view of MolTrack EMS", () => openEmsView(page, "Compounds"));
    enter(page, "MolTrack EMS");
    await session.step(16, "When user clicks on more actions in app ribbon", () => clickOn(page, el("more actions in app ribbon")), undefined, "Bulk registration is the Import of the Compounds view");
    await session.step(17, "And user clicks on \"Import…\" menu item", () => clickOn(page, el("\"Import…\" menu item")));
    await session.step(18, "Then \"Register compounds\" dialog should be visible", () => shouldBe(page, el("\"Register compounds\" dialog"), "visible"));
    await session.step(20, "And \"Data\" input in \"Register compounds\" dialog should contain text \"bulk-registration\"", () => shouldContainText(page, el("\"Data\" input in \"Register compounds\" dialog"), "bulk-registration"), undefined, "The table to register");
    await session.step(22, "And \"Register compounds\" dialog should contain text \"5 rows: 3 new, 1 already registered, 1 cannot be read\"", () => shouldContainText(page, el("\"Register compounds\" dialog"), "5 rows: 3 new, 1 already registered, 1 cannot be read"), undefined, "The preview: each row's structure and what it becomes");
    await session.step(23, "When user clicks on \"OK\" button in \"Register compounds\" dialog", () => clickOn(page, el("\"OK\" button in \"Register compounds\" dialog")));
    await session.step(25, "Then \"Import compounds\" dialog should be visible", () => shouldBe(page, el("\"Import compounds\" dialog"), "visible"), undefined, "The new compounds go on to the import: columns, preview, report");
    await session.step(26, "When user clicks on \"Next\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Next\" button in \"Import compounds\" dialog")));
    await session.step(28, "Then import mapping should contain text \"Common Name\"", () => shouldContainText(page, el("import mapping"), "Common Name"), undefined, "Each column of the table, and the compound column it fills");
    await session.step(29, "When user clicks on \"Next\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Next\" button in \"Import compounds\" dialog")));
    await session.step(31, "Then import preview should contain text \"3 will be added\"", () => shouldContainText(page, el("import preview"), "3 will be added"), undefined, "The preview: three new compounds");
    await session.step(32, "When user clicks on \"Next\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Next\" button in \"Import compounds\" dialog")));
    await session.step(34, "Then import report should contain text \"3 inserted\"", () => shouldContainText(page, el("import report"), "3 inserted"), undefined, "The report");
    await session.step(35, "When user clicks on \"Close\" button in \"Import compounds\" dialog", () => clickOn(page, el("\"Close\" button in \"Import compounds\" dialog")));
    await session.step(37, "Then \"Grok ID\" column of table \"bulk-registration\" should have type \"string\"", () => columnTypeOfTable(page, "Grok ID", "bulk-registration", "string"), undefined, "Every row of the table gets its Grok ID back");
    await session.step(38, "And the bulk registration example should have 3 new compounds on the server", () => bulkRegistered(page, 3));
  });
});
