/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/open-a-jmp-file.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableViewOpen, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Open a JMP file in the browser", () => {
  const session = feature(test, "features/guides/open-a-jmp-file.feature", import.meta.url);
  test("Double-click a JMP file to open it with its saved analyses", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And no table view is open", () => noTableViewOpen(page));
    await session.step(13, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(14, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(15, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(16, "And Files---App-Data---JmpReader tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader tree node inside browse tree")));
    await session.step(17, "And Files---App-Data---JmpReader---public tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---JmpReader---public tree node inside browse tree")));
    await session.step(18, "When user double-clicks on Files---App-Data---JmpReader---public---mdpi-bioengineering-11-01089-wp1.jmp tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---JmpReader---public---mdpi-bioengineering-11-01089-wp1.jmp tree node inside browse tree")));
    await session.step(19, "Then the \"mdpi-bioengineering-11-01089-wp1\" table view should open with 21 rows", () => tableViewOpened(page, "mdpi-bioengineering-11-01089-wp1", 21));
    await session.step(20, "And the open tableview should have 3 scatter plot viewers", () => viewerCount(page, 3, "scatter plot"));
    await session.step(21, "When user clicks on \"About this JMP file\" icon", () => clickOn(page, el("\"About this JMP file\" icon")));
    await session.step(22, "Then \"Saved analyses\" section should be visible", () => shouldBe(page, el("\"Saved analyses\" section"), "visible"));
    await session.step(23, "And \"Saved analyses\" section should contain text \"Scatter plot\"", () => shouldContainText(page, el("\"Saved analyses\" section"), "Scatter plot"));
    await session.step(24, "When user clicks on CANCEL button", () => clickOn(page, el("CANCEL button")));
    await session.step(25, "Then the \"mdpi-bioengineering-11-01089-wp1.jmp\" dialog should close", () => dialogCloses(page, "mdpi-bioengineering-11-01089-wp1.jmp"));
  });
});
