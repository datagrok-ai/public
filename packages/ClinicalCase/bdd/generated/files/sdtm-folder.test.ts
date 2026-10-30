/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/files/sdtm-folder.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, currentViewShows} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A folder of SDTM data in Browse > Files offers to run Clinical Case", () => {
  const session = feature(test, "features/files/sdtm-folder.feature", import.meta.url);
  test("The study folder shows the launcher with the study id from its dm.csv", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(11, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(12, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(13, "And Files---App-Data---ClinicalCase tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---ClinicalCase tree node inside browse tree")));
    await session.step(14, "And Files---App-Data---ClinicalCase---SDTM tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---ClinicalCase---SDTM tree node inside browse tree")));
    await session.step(17, "When user clicks on Files---App-Data---ClinicalCase---SDTM---CDISCPILOT01 tree node inside browse tree", () => clickOn(page, el("Files---App-Data---ClinicalCase---SDTM---CDISCPILOT01 tree node inside browse tree")));
    await session.step(18, "Then the current view should show the text \"Folder contains SDTM data\"", () => currentViewShows(page, "Folder contains SDTM data"));
    await session.step(19, "And the current view should show the text \"Study ID: CDISCPILOT01\"", () => currentViewShows(page, "Study ID: CDISCPILOT01"));
    await session.step(20, "And \"Run ClinicalCase\" button should be visible", () => shouldBe(page, el("\"Run ClinicalCase\" button"), "visible"));
    await session.step(21, "And no errors should have been logged", () => noErrors(page));
  });
});
