/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/app/clinical-case-app.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, isExpanded, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeCurrentView, currentViewShows, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Clinical Case opens from Browse > Apps with its demo studies and the Import study view", () => {
  const session = feature(test, "features/app/clinical-case-app.feature", import.meta.url);
  test("Clinical Case opens from Browse > Apps with its demo studies and the Import study view", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await run.scenario("The start page lists the demo studies with their subject counts", async () => {
      await session.step(21, "When user clicks on Apps---Clinical-Case tree node inside browse tree", () => clickOn(page, el("Apps---Clinical-Case tree node inside browse tree")));
      await session.step(22, "Then the \"Clinical Case\" view should be current", () => viewIsCurrent(page, "Clinical Case"));
      await session.step(23, "And the current view should show the text \"Visualize and explore your SDTM data\"", () => currentViewShows(page, "Visualize and explore your SDTM data"));
      await session.step(24, "And the current view should show the text \"Study Standard Protocol Total subjects Start date End date CDISCPILOT01 SDTM 306 nida-ctn-0001 SDTM 138\"", () => currentViewShows(page, "Study Standard Protocol Total subjects Start date End date CDISCPILOT01 SDTM 306 nida-ctn-0001 SDTM 138"));
      await session.step(25, "And \"Import study...\" button should be visible", () => shouldBe(page, el("\"Import study...\" button"), "visible"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Clinical Case node lists the studies and Import study", async () => {
      await session.step(29, "When user expands Apps---Clinical-Case tree node inside browse tree", () => expand(page, el("Apps---Clinical-Case tree node inside browse tree")));
      await session.step(30, "Then Apps---Clinical-Case---CDISCPILOT01 tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Clinical-Case---CDISCPILOT01 tree node inside browse tree"), "visible"));
      await session.step(31, "And Apps---Clinical-Case---nida-ctn-0001 tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Clinical-Case---nida-ctn-0001 tree node inside browse tree"), "visible"));
      await session.step(32, "And Apps---Clinical-Case---Import-study tree node inside browse tree should be visible", () => shouldBe(page, el("Apps---Clinical-Case---Import-study tree node inside browse tree"), "visible"));
      await session.step(33, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Import study... opens the import view, which refuses files without a study definition", async () => {
      await session.step(36, "When user clicks on Apps---Clinical-Case tree node inside browse tree", () => clickOn(page, el("Apps---Clinical-Case tree node inside browse tree")));
      await session.step(37, "Then the \"Clinical Case\" view should be current", () => viewIsCurrent(page, "Clinical Case"));
      await session.step(38, "When user clicks on \"Import study...\" button", () => clickOn(page, el("\"Import study...\" button")));
      await session.step(39, "Then the \"Import Study - Clinical Case\" view should be current", () => viewIsCurrent(page, "Import Study - Clinical Case"));
      await session.step(40, "And the current view should show the text \"Import study files\"", () => currentViewShows(page, "Import study files"));
      await session.step(41, "And \"Import\" button should be disabled", () => shouldBe(page, el("\"Import\" button"), "disabled"));
      await session.step(42, "When user uploads \"../files/SDTM/CDISCPILOT01/dm.csv\" through \"Open files\" icon", () => uploadThrough(page, "../files/SDTM/CDISCPILOT01/dm.csv", el("\"Open files\" icon")));
      await session.step(43, "Then the current view should show the text \"Error: Invalid or missing define.xml/study.json\"", () => currentViewShows(page, "Error: Invalid or missing define.xml/study.json"));
      await session.step(44, "And the current view should show the text \"dm.csv\"", () => currentViewShows(page, "dm.csv"));
      await session.step(45, "And an error balloon containing \"Invalid or missing define.xml/study.json\" should have been shown", () => errorBalloonText(page, "Invalid or missing define.xml/study.json"));
      await session.step(46, "And \"Import\" button should be disabled", () => shouldBe(page, el("\"Import\" button"), "disabled"));
      await session.step(47, "When user closes the current view", () => closeCurrentView(page));
      await session.step(48, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
