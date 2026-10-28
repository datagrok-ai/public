/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/preview.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [jmp.preview]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {above} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, currentViewShows, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableViewOpen, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A JMP file in the Browse tree previews before it is opened", () => {
  const session = feature(test, "features/preview.feature", import.meta.url);
  test("The preview tells what the file is, and Open opens the table", {tag: ["@realizes:jmp.preview"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(17, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(18, "And Files---App-Data---JmpReader tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader tree node inside browse tree")));
    await session.step(21, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(22, "And Files---App-Data---JmpReader---jmpreader tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader---jmpreader tree node inside browse tree")));
    await session.step(23, "When user double-clicks on Files---App-Data---JmpReader---jmpreader---example1.jmp tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---JmpReader---jmpreader---example1.jmp tree node inside browse tree")));
    await session.step(24, "Then the \"example1.jmp\" view should be current", () => viewIsCurrent(page, "example1.jmp"));
    await session.step(25, "And the current view should show the text \"JMP 17.2.0\"", () => currentViewShows(page, "JMP 17.2.0"));
    await session.step(26, "And the current view should show the text \"4 rows × 12 columns\"", () => currentViewShows(page, "4 rows × 12 columns"));
    await session.step(27, "And \"Columns not imported\" section should be absent", () => shouldBe(page, el("\"Columns not imported\" section"), "absent"));
    await session.step(28, "When user clicks on Open button", () => clickOn(page, el("Open button")));
    await session.step(29, "Then the \"example1\" table view should open with 4 rows", () => tableViewOpened(page, "example1", 4));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The preview lists the columns it can't import", {tag: ["@realizes:jmp.preview"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(17, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(18, "And Files---App-Data---JmpReader tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader tree node inside browse tree")));
    await session.step(34, "Given Files---App-Data---JmpReader---tests tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader---tests tree node inside browse tree")));
    await session.step(35, "When user double-clicks on Files---App-Data---JmpReader---tests---unsupported-column.jmp tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---JmpReader---tests---unsupported-column.jmp tree node inside browse tree")));
    await session.step(36, "Then the \"unsupported-column.jmp\" view should be current", () => viewIsCurrent(page, "unsupported-column.jmp"));
    await session.step(37, "And the current view should show the text \"4 rows × 11 columns\"", () => currentViewShows(page, "4 rows × 11 columns"));
    await session.step(38, "And \"Columns not imported\" section should be visible", () => shouldBe(page, el("\"Columns not imported\" section"), "visible"));
    await session.step(39, "And \"Columns not imported\" section should contain text \"1 of 12\"", () => shouldContainText(page, el("\"Columns not imported\" section"), "1 of 12"));
    await session.step(40, "And \"Columns not imported\" section should be above grid", () => above(page, el("\"Columns not imported\" section"), el("grid")));
    await session.step(41, "When user expands \"Columns not imported\" section", () => expand(page, el("\"Columns not imported\" section")));
    await session.step(42, "Then \"Columns not imported\" section should contain text \"formula\"", () => shouldContainText(page, el("\"Columns not imported\" section"), "formula"));
    await session.step(43, "And \"Columns not imported\" section should contain text \"Unsupported column layout\"", () => shouldContainText(page, el("\"Columns not imported\" section"), "Unsupported column layout"));
    await session.step(44, "And \"Columns not imported\" section should be above grid", () => above(page, el("\"Columns not imported\" section"), el("grid")));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
