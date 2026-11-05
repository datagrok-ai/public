/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/open.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [jmp.import]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, menuFileChooser, noTableViewOpen, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A JMP data table from disk opens as a table", () => {
  const session = feature(test, "features/open.feature", import.meta.url);
  test("The table opens when dropped onto the window", {tag: ["@realizes:jmp.import"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "When user drops the \"fixtures/example1.jmp\" file of the project onto browse tree", () => dropFile(page, "fixtures/example1.jmp", el("browse tree")));
    await session.step(24, "Then the \"example1\" table view should open with 4 rows", () => tableViewOpened(page, "example1", 4));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
    await session.step(26, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The table opens through the Open local file icon and the file chooser", {tag: ["@realizes:jmp.import"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "When user uploads \"fixtures/example1.jmp\" through \"Open local file\" icon in browse toolbar", () => uploadThrough(page, "fixtures/example1.jmp", el("\"Open local file\" icon in browse toolbar")));
    await session.step(32, "Then the \"example1\" table view should open with 4 rows", () => tableViewOpened(page, "example1", 4));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
    await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The table opens through File > Open > File... and the file chooser", {tag: ["@realizes:jmp.import"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(37, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(38, "When user clicks on Menu tab", () => clickOn(page, el("Menu tab")));
    await session.step(39, "And user picks \"File > Open > File...\" from the open menu and chooses the \"fixtures/example1.jmp\" file of the project", () => menuFileChooser(page, "File > Open > File...", "fixtures/example1.jmp"));
    await session.step(40, "Then the \"example1\" table view should open with 4 rows", () => tableViewOpened(page, "example1", 4));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
