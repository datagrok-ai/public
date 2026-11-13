/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/plates.feature
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
import {clearField, clickOn, enterInto, isExpanded, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelOpen, packageInstalled, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Apps > Plates: the search form follows its template, a template opens in its editor, Create warns with no plate", () => {
  const session = feature(test, "features/apps/plates.feature", import.meta.url);
  test("Apps > Plates: the search form follows its template, a template opens in its editor, Create warns with no plate", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the \"Plates\" package is installed", () => packageInstalled(page, "Plates"));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(18, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
    await run.scenario("Apps > Plates opens the list of plates with their properties as columns", async () => {
      await session.step(21, "When user clicks on Apps---Plates tree node inside browse tree", () => clickOn(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(22, "Then the \"Plates\" view should be current", () => viewIsCurrent(page, "Plates"));
      await session.step(23, "And the table should have a column \"Imaging device\"", () => hasColumn(page, "Imaging device"));
      await session.step(24, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Search Wells checks a numeric criterion and swaps its well inputs for the template picked", async () => {
      await session.step(27, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(28, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(29, "When user clicks on Apps---Plates---Search-wells tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Search-wells tree node inside browse tree")));
      await session.step(30, "Then the \"Search Wells\" view should be current", () => viewIsCurrent(page, "Search Wells"));
      await session.step(31, "And \"Well cell count\" input should be visible", () => shouldBe(page, el("\"Well cell count\" input"), "visible"));
      await session.step(32, "When user enters \"abc\" into \"Well cell count\" input", () => enterInto(page, "abc", el("\"Well cell count\" input")));
      await session.step(33, "Then \"Well cell count\" input should be invalid", () => shouldBe(page, el("\"Well cell count\" input"), "invalid"));
      await session.step(34, "When user clears \"Well cell count\" input", () => clearField(page, el("\"Well cell count\" input")));
      await session.step(35, "Then \"Well cell count\" input should be valid", () => shouldBe(page, el("\"Well cell count\" input"), "valid"));
      await session.step(36, "When user selects \"Dose-response\" in \"Template\" input", () => selectIn(page, "Dose-response", el("\"Template\" input")));
      await session.step(37, "Then \"Concentration\" input should be visible", () => shouldBe(page, el("\"Concentration\" input"), "visible"));
      await session.step(38, "And \"Well cell count\" input should be hidden", () => shouldBe(page, el("\"Well cell count\" input"), "hidden"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("EDIT opens the template in its editor", async () => {
      await session.step(42, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(43, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(44, "When user clicks on Apps---Plates---Templates tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Templates tree node inside browse tree")));
      await session.step(45, "Then the \"Templates\" view should be current", () => viewIsCurrent(page, "Templates"));
      await session.step(46, "When user clicks on first \"EDIT\" button", () => clickOn(page, el("first \"EDIT\" button")));
      await session.step(47, "Then the \"Templates / Cell counting\" view should be current", () => viewIsCurrent(page, "Templates / Cell counting"));
      await session.step(48, "And \"Name\" input should have value \"Cell counting\"", () => shouldHaveValue(page, el("\"Name\" input"), "Cell counting"));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("CREATE with no plate imported only warns", async () => {
      await session.step(52, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(53, "And Apps---Plates tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Plates tree node inside browse tree")));
      await session.step(54, "When user clicks on Apps---Plates---Create tree node inside browse tree", () => clickOn(page, el("Apps---Plates---Create tree node inside browse tree")));
      await session.step(55, "Then the \"Create Plate\" view should be current", () => viewIsCurrent(page, "Create Plate"));
      await session.step(56, "When user clicks on \"CREATE\" button", () => clickOn(page, el("\"CREATE\" button")));
      await session.step(57, "Then a warning balloon containing \"No active plate is selected.\" should have been shown", () => warningBalloonText(page, "No active plate is selected."));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
