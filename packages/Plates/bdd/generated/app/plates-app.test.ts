/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/app/plates-app.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, isExpanded, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewShows, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingIs, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Apps > Plates: the plate list, the searches, Templates and Create Plate open with their content", () => {
  const session = feature(test, "features/app/plates-app.feature", import.meta.url);
  test("Apps > Plates: the plate list, the searches, Templates and Create Plate open with their content", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(26, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
    await run.scenario("Apps > Plates lists the registered plates with their plate properties as columns", async () => {
      await session.step(29, "When user clicks on Apps---Plates tree node inside browse tree", () => clickOn(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(30, "Then the \"Plates\" view should be current", () => viewIsCurrent(page, "Plates"));
      await session.step(31, "And the table should have 40 rows", () => rowCount(page, 40));
      await session.step(32, "And \"barcode\" column should have semantic type \"Plate Barcode\"", () => columnSemType(page, "barcode", "Plate Barcode"));
      await session.step(33, "And the table should have a column \"Imaging device\"", () => hasColumn(page, "Imaging device"));
      await session.step(34, "And the table should have a column \"Project\"", () => hasColumn(page, "Project"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Search Wells offers the well properties of the template picked in it", async () => {
      await session.step(38, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(39, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(40, "When user clicks on Apps---Plates---Search-wells tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Search-wells tree node inside browse tree")));
      await session.step(41, "Then the \"Search Wells\" view should be current", () => viewIsCurrent(page, "Search Wells"));
      await session.step(42, "And \"Template\" input should have value \"Cell counting\"", () => shouldHaveValue(page, el("\"Template\" input"), "Cell counting"));
      await session.step(43, "And \"Well cell count\" input should be visible", () => shouldBe(page, el("\"Well cell count\" input"), "visible"));
      await session.step(44, "When user selects \"Dose-response\" in \"Template\" input", () => selectIn(page, "Dose-response", el("\"Template\" input")));
      await session.step(45, "Then \"Concentration\" input should be visible", () => shouldBe(page, el("\"Concentration\" input"), "visible"));
      await session.step(46, "And \"Activity\" input should be visible", () => shouldBe(page, el("\"Activity\" input"), "visible"));
      await session.step(47, "And \"Well cell count\" input should be hidden", () => shouldBe(page, el("\"Well cell count\" input"), "hidden"));
      await session.step(48, "When user selects \"Cell counting\" in \"Template\" input", () => selectIn(page, "Cell counting", el("\"Template\" input")));
      await session.step(49, "Then \"Well cell count\" input should be visible", () => shouldBe(page, el("\"Well cell count\" input"), "visible"));
      await session.step(50, "And \"Concentration\" input should be hidden", () => shouldBe(page, el("\"Concentration\" input"), "hidden"));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Search Plates narrows the plates found by a number typed into a plate property", async () => {
      await session.step(54, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(55, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(56, "When user clicks on Apps---Plates---Search-plates tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Search-plates tree node inside browse tree")));
      await session.step(57, "Then the \"Search Plates\" view should be current", () => viewIsCurrent(page, "Search Plates"));
      await session.step(58, "And the \"rows\" reading of grid should be 40", () => readingIs(page, "rows", el("grid"), 40));
      await session.step(59, "When user enters \">=0\" into \"Plate cell count\" input", () => enterInto(page, ">=0", el("\"Plate cell count\" input")));
      await session.step(60, "Then the \"rows\" reading of grid should be 20", () => readingIs(page, "rows", el("grid"), 20));
      await session.step(61, "When user clears \"Plate cell count\" input", () => clearField(page, el("\"Plate cell count\" input")));
      await session.step(62, "Then the \"rows\" reading of grid should be 40", () => readingIs(page, "rows", el("grid"), 40));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Templates lists the templates; EDIT opens one, CREATE NEW a new one to save", async () => {
      await session.step(66, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(67, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(68, "When user clicks on Apps---Plates---Templates tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Templates tree node inside browse tree")));
      await session.step(69, "Then the \"Templates\" view should be current", () => viewIsCurrent(page, "Templates"));
      await session.step(70, "And the current view should show the text \"Cell counting\"", () => currentViewShows(page, "Cell counting"));
      await session.step(71, "And the current view should show the text \"Dose-response\"", () => currentViewShows(page, "Dose-response"));
      await session.step(72, "When user clicks on first \"EDIT\" button", () => clickOn(page, el("first \"EDIT\" button")));
      await session.step(73, "Then the \"Templates / Cell counting\" view should be current", () => viewIsCurrent(page, "Templates / Cell counting"));
      await session.step(74, "And \"Name\" input should have value \"Cell counting\"", () => shouldHaveValue(page, el("\"Name\" input"), "Cell counting"));
      await session.step(75, "And \"SAVE\" button should be hidden", () => shouldBe(page, el("\"SAVE\" button"), "hidden"));
      await session.step(76, "When user clicks on Apps---Plates---Search-wells tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Search-wells tree node inside browse tree")));
      await session.step(77, "And user clicks on Apps---Plates---Templates tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Templates tree node inside browse tree")));
      await session.step(78, "Then the \"Templates\" view should be current", () => viewIsCurrent(page, "Templates"));
      await session.step(79, "When user clicks on \"CREATE NEW\" button", () => clickOn(page, el("\"CREATE NEW\" button")));
      await session.step(80, "Then the \"Templates / New Template\" view should be current", () => viewIsCurrent(page, "Templates / New Template"));
      await session.step(81, "And \"Name\" input should have value \"New Template\"", () => shouldHaveValue(page, el("\"Name\" input"), "New Template"));
      await session.step(82, "And \"SAVE\" button should be visible", () => shouldBe(page, el("\"SAVE\" button"), "visible"));
      await session.step(83, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Create Plate shows an empty plate with an analysis tab each, and CREATE with no plate imported only warns", async () => {
      await session.step(86, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(87, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(88, "When user clicks on Apps---Plates---Create tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Create tree node inside browse tree")));
      await session.step(89, "Then the \"Create Plate\" view should be current", () => viewIsCurrent(page, "Create Plate"));
      await session.step(90, "And \"Plate View\" tab should be visible", () => shouldBe(page, el("\"Plate View\" tab"), "visible"));
      await session.step(91, "And \"qPCR Analysis\" tab should be visible", () => shouldBe(page, el("\"qPCR Analysis\" tab"), "visible"));
      await session.step(92, "When user clicks on \"Dose Response\" tab", () => clickOn(page, el("\"Dose Response\" tab")));
      await session.step(93, "Then the current view should show the text \"To see dose response curves, please map the following required properties: Activity, Concentration, SampleID.\"", () => currentViewShows(page, "To see dose response curves, please map the following required properties: Activity, Concentration, SampleID."));
      await session.step(94, "When user clicks on \"CREATE\" button", () => clickOn(page, el("\"CREATE\" button")));
      await session.step(95, "Then a warning balloon containing \"No active plate is selected.\" should have been shown", () => warningBalloonText(page, "No active plate is selected."));
      await session.step(96, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
