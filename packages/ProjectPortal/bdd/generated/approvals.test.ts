/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/approvals.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProjectAtStage} from '../bindings/steps.js';
import {clickOn, enterInto, pressKey, shouldBe, shouldContainText, shouldNotBe, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Approval of gate decisions", () => {
  const session = feature(test, "features/approvals.feature", import.meta.url);
  test("An approver proposes an advance, which waits for another approver", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given a scratch project \"BDD-APPROVAL\" at stage \"TI\" in the \"EX-BIOD\" portfolio", () => scratchProjectAtStage(page, "BDD-APPROVAL", "TI", "EX-BIOD"));
    await session.step(12, "When user clicks on biologics discovery portfolio", () => clickOn(page, el("biologics discovery portfolio")));
    await session.step(13, "And user clicks on approval project card", () => clickOn(page, el("approval project card")));
    await session.step(14, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(15, "Then approval setting should contain text \"Needs 1 approver besides the proposer: advance, terminate\"", () => shouldContainText(page, el("approval setting"), "Needs 1 approver besides the proposer: advance, terminate"));
    await session.step(16, "When user clicks on binding criterion met button", () => clickOn(page, el("binding criterion met button")));
    await session.step(17, "Then advance button should contain text \"Propose advance to Lead Generation\"", () => shouldContainText(page, el("advance button"), "Propose advance to Lead Generation"));
    await session.step(18, "When user enters \"Binders below 10 nM\" into rationale input", () => enterInto(page, "Binders below 10 nM", el("rationale input")));
    await session.step(19, "And user clicks on advance button", () => clickOn(page, el("advance button")));
    await session.step(20, "Then pending approvals should contain text \"BDD-APPROVAL: advance to Lead Generation\"", () => shouldContainText(page, el("pending approvals"), "BDD-APPROVAL: advance to Lead Generation"));
    await session.step(21, "And pending approvals should contain text \"0 of 1 approval\"", () => shouldContainText(page, el("pending approvals"), "0 of 1 approval"));
    await session.step(22, "And approval result should contain text \"proposed; waits for 1 approval\"", () => shouldContainText(page, el("approval result"), "proposed; waits for 1 approval"));
    await session.step(23, "And advance button should be disabled", () => shouldBe(page, el("advance button"), "disabled"));
    await session.step(24, "And first approval approve button should not be visible", () => shouldNotBe(page, el("first approval approve button"), "visible"));
    await session.step(25, "And first approval withdraw button should be visible", () => shouldBe(page, el("first approval withdraw button"), "visible"));
    await session.step(26, "When user clicks on waiting for me button", () => clickOn(page, el("waiting for me button")));
    await session.step(27, "Then waiting dialog should be visible", () => shouldBe(page, el("waiting dialog"), "visible"));
    await session.step(28, "And waiting list should not contain text \"BDD-APPROVAL\"", () => shouldNotContainText(page, el("waiting list"), "BDD-APPROVAL"));
    await session.step(29, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(30, "And user clicks on first approval withdraw button", () => clickOn(page, el("first approval withdraw button")));
    await session.step(31, "Then no pending approvals should be visible", () => shouldBe(page, el("no pending approvals"), "visible"));
    await session.step(32, "And closed approvals should contain text \"withdrawn\"", () => shouldContainText(page, el("closed approvals"), "withdrawn"));
    await session.step(33, "And advance button should be enabled", () => shouldBe(page, el("advance button"), "enabled"));
  });
});
