/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/problem-files.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [jmp.import]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {columnCount, hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelOpen, currentViewShows, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, noTableLeft, noTableViewOpen, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A JMP file the reader can't fully read says so", () => {
  const session = feature(test, "features/problem-files.feature", import.meta.url);
  test("An unreadable column is left out with a warning, the rest opens", {tag: ["@realizes:jmp.import"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "When user drops the \"fixtures/unsupported-column.jmp\" file of the project onto browse tree", () => dropFile(page, "fixtures/unsupported-column.jmp", el("browse tree")));
    await session.step(18, "Then the \"unsupported-column\" table view should open with 4 rows", () => tableViewOpened(page, "unsupported-column", 4));
    await session.step(19, "And the table should have 11 columns", () => columnCount(page, 11));
    await session.step(20, "And the table should not have a column \"formula\"", () => hasNoColumn(page, "formula"));
    await session.step(21, "And a warning balloon containing \"1 column not imported: formula\" should have been shown", () => warningBalloonText(page, "1 column not imported: formula"));
    await session.step(22, "And no errors should have been logged", () => noErrors(page));
  });
  test("A file that isn't a JMP table says why, and opens no table", {tag: ["@realizes:jmp.import"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "When user drops the \"fixtures/not-a-jmp.jmp\" file of the project onto browse tree", () => dropFile(page, "fixtures/not-a-jmp.jmp", el("browse tree")));
    await session.step(28, "Then the \"not-a-jmp.jmp\" view should be current", () => viewIsCurrent(page, "not-a-jmp.jmp"));
    await session.step(29, "And the current view should show the text \"This file can't be opened as a JMP table\"", () => currentViewShows(page, "This file can't be opened as a JMP table"));
    await session.step(30, "And the current view should show the text \"Not a JMP data table\"", () => currentViewShows(page, "Not a JMP data table"));
    await session.step(31, "And no table should be left in the workspace", () => noTableLeft(page));
    await session.step(32, "And no errors should have been logged", () => noErrors(page));
  });
});
