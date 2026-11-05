/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/files/plate-text-file.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, closeCurrentView, refreshBrowse, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInUsersFiles} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A plate reader's text file opens as a plate, from the computer and from Browse > Files", () => {
  const session = feature(test, "features/files/plate-text-file.feature", import.meta.url);
  test("A plate reader's text file opens as a plate, from the computer and from Browse > Files", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("A BMG PHERAstar file opened from the computer comes up as a plate with a tab per channel", async () => {
      await session.step(24, "When user uploads \"fixtures/bdd-bmg-plate.txt\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/bdd-bmg-plate.txt", el("\"Open local file\" icon inside browse toolbar")));
      await session.step(25, "Then the \"bdd-bmg-plate.txt\" view should be current", () => viewIsCurrent(page, "bdd-bmg-plate.txt"));
      await session.step(26, "And \"Summary\" tab should be visible", () => shouldBe(page, el("\"Summary\" tab"), "visible"));
      await session.step(27, "And \"Channel A\" tab should be visible", () => shouldBe(page, el("\"Channel A\" tab"), "visible"));
      await session.step(28, "And \"Channel B\" tab should be visible", () => shouldBe(page, el("\"Channel B\" tab"), "visible"));
      await session.step(29, "And the \"rows\" reading of heat map viewer should be 8", () => readingIs(page, "rows", el("heat map viewer"), 8));
      await session.step(30, "When user clicks on \"Channel A\" tab", () => clickOn(page, el("\"Channel A\" tab")));
      await session.step(31, "Then the \"text of cell 1 of 1\" reading of grid should be \"1000.00\"", () => readingReads(page, "text of cell 1 of 1", el("grid"), "1000.00"));
      await session.step(32, "And the \"text of cell 2 of 1\" reading of grid should be \"1200.00\"", () => readingReads(page, "text of cell 2 of 1", el("grid"), "1200.00"));
      await session.step(33, "And the \"text of cell 1 of 2\" reading of grid should be \"1100.00\"", () => readingReads(page, "text of cell 1 of 2", el("grid"), "1100.00"));
      await session.step(34, "And the \"text of cell 8 of 12\" reading of grid should be \"1300.00\"", () => readingReads(page, "text of cell 8 of 12", el("grid"), "1300.00"));
      await session.step(35, "And the \"text of cell 2 of 2\" reading of grid should be \"\"", () => readingReads(page, "text of cell 2 of 2", el("grid"), ""));
      await session.step(36, "When user clicks on \"Channel B\" tab", () => clickOn(page, el("\"Channel B\" tab")));
      await session.step(37, "Then the \"text of cell 1 of 1\" reading of grid should be \"2000.00\"", () => readingReads(page, "text of cell 1 of 1", el("grid"), "2000.00"));
      await session.step(38, "And the \"text of cell 8 of 12\" reading of grid should be \"2300.00\"", () => readingReads(page, "text of cell 8 of 12", el("grid"), "2300.00"));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
      await session.step(40, "When user closes the current view", () => closeCurrentView(page));
      await session.step(41, "Then \"bdd-bmg-plate.txt\" view tab should be absent", () => shouldBe(page, el("\"bdd-bmg-plate.txt\" view tab"), "absent"));
    });
    await run.scenario("The same file in Browse > Files previews as a plate", async () => {
      await session.step(44, "Given the \"fixtures/bdd-bmg-plate.txt\" file of the project is in the user's files as \"bdd-plate-{time}.txt\"", () => fixtureInUsersFiles(page, "fixtures/bdd-bmg-plate.txt", session.text("bdd-plate-{time}.txt")));
      await session.step(45, "When user closes all views", () => closeAllViews(page));
      await session.step(46, "And user refreshes the browse tree", () => refreshBrowse(page));
      await session.step(47, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(48, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
      await session.step(49, "And user clicks on Files---My-files---bdd-plate-{time}.txt tree node inside browse tree", () => clickOn(page, el(session.text("Files---My-files---bdd-plate-{time}.txt tree node inside browse tree"))));
      await session.step(50, "Then the \"bdd-plate-{time}.txt\" view should be current", () => viewIsCurrent(page, session.text("bdd-plate-{time}.txt")));
      await session.step(51, "And the \"rows\" reading of heat map viewer should be 8", () => readingIs(page, "rows", el("heat map viewer"), 8));
      await session.step(52, "When user clicks on \"Channel A\" tab", () => clickOn(page, el("\"Channel A\" tab")));
      await session.step(53, "Then the \"text of cell 1 of 1\" reading of grid should be \"1000.00\"", () => readingReads(page, "text of cell 1 of 1", el("grid"), "1000.00"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
