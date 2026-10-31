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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableViewOpen, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A JMP file in the Browse tree opens as its table", () => {
  const session = feature(test, "features/preview.feature", import.meta.url);
  test("Double-click opens the table, and the summary tells what the file is", {tag: ["@realizes:jmp.preview"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(17, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(18, "And Files---App-Data---JmpReader tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader tree node inside browse tree")));
    await session.step(21, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(22, "And Files---App-Data---JmpReader---jmpreader tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader---jmpreader tree node inside browse tree")));
    await session.step(23, "When user double-clicks on Files---App-Data---JmpReader---jmpreader---example1.jmp tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---JmpReader---jmpreader---example1.jmp tree node inside browse tree")));
    await session.step(24, "Then the \"example1\" table view should open with 4 rows", () => tableViewOpened(page, "example1", 4));
    await session.step(25, "When user clicks on \"About this JMP file\" icon", () => clickOn(page, el("\"About this JMP file\" icon")));
    await session.step(26, "Then \"example1.jmp\" dialog should contain text \"JMP 17.2.0\"", () => shouldContainText(page, el("\"example1.jmp\" dialog"), "JMP 17.2.0"));
    await session.step(27, "And \"example1.jmp\" dialog should contain text \"4 rows × 12 columns\"", () => shouldContainText(page, el("\"example1.jmp\" dialog"), "4 rows × 12 columns"));
    await session.step(28, "And \"Columns not imported\" section should be absent", () => shouldBe(page, el("\"Columns not imported\" section"), "absent"));
    await session.step(29, "When user clicks on CANCEL button", () => clickOn(page, el("CANCEL button")));
    await session.step(30, "Then the \"example1.jmp\" dialog should close", () => dialogCloses(page, "example1.jmp"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
    await session.step(32, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The columns it can't import are named, and listed in the summary", {tag: ["@realizes:jmp.preview"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(17, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(18, "And Files---App-Data---JmpReader tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader tree node inside browse tree")));
    await session.step(35, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(36, "And Files---App-Data---JmpReader---tests tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader---tests tree node inside browse tree")));
    await session.step(37, "When user double-clicks on Files---App-Data---JmpReader---tests---unsupported-column.jmp tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---JmpReader---tests---unsupported-column.jmp tree node inside browse tree")));
    await session.step(38, "Then the \"unsupported-column\" table view should open with 4 rows", () => tableViewOpened(page, "unsupported-column", 4));
    await session.step(39, "And a warning balloon containing \"1 column not imported: formula\" should have been shown", () => warningBalloonText(page, "1 column not imported: formula"));
    await session.step(40, "When user clicks on \"About this JMP file\" icon", () => clickOn(page, el("\"About this JMP file\" icon")));
    await session.step(41, "Then \"unsupported-column.jmp\" dialog should contain text \"4 rows × 11 columns\"", () => shouldContainText(page, el("\"unsupported-column.jmp\" dialog"), "4 rows × 11 columns"));
    await session.step(42, "And \"Columns not imported\" section should be visible", () => shouldBe(page, el("\"Columns not imported\" section"), "visible"));
    await session.step(43, "And \"Columns not imported\" section should contain text \"1 of 12\"", () => shouldContainText(page, el("\"Columns not imported\" section"), "1 of 12"));
    await session.step(44, "When user expands \"Columns not imported\" section", () => expand(page, el("\"Columns not imported\" section")));
    await session.step(45, "Then \"Columns not imported\" section should contain text \"formula\"", () => shouldContainText(page, el("\"Columns not imported\" section"), "formula"));
    await session.step(46, "And \"Columns not imported\" section should contain text \"Unsupported column layout\"", () => shouldContainText(page, el("\"Columns not imported\" section"), "Unsupported column layout"));
    await session.step(47, "When user clicks on CANCEL button", () => clickOn(page, el("CANCEL button")));
    await session.step(48, "Then the \"unsupported-column.jmp\" dialog should close", () => dialogCloses(page, "unsupported-column.jmp"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
