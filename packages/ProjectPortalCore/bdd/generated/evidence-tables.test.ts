/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/evidence-tables.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProjectAtStage, scratchWorkspace, titledRunSaved} from '../bindings/steps.js';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldNotBe, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Evidence tables generated into documents", () => {
  const session = feature(test, "features/evidence-tables.feature", import.meta.url);
  test("The process description's evidence is generated into a draft document", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given a scratch project \"BDD-ET\" at stage \"D4\" in the \"EX-BIOC\" portfolio", () => scratchProjectAtStage(page, "BDD-ET", "D4", "EX-BIOC"));
    await session.step(13, "And the workspace of the scratch project \"BDD-ET\" in the \"EX-BIOC\" portfolio is built", () => scratchWorkspace(page, "BDD-ET", "EX-BIOC"));
    await session.step(14, "And a run of the example model at feed rate 2 titled \"Feed 2\" saved to the scratch project \"BDD-ET\" in the \"EX-BIOC\" portfolio", () => titledRunSaved(page, 2, "Feed 2", "BDD-ET", "EX-BIOC"));
    await session.step(15, "And a run of the example model at feed rate 4 titled \"Feed 4\" saved to the scratch project \"BDD-ET\" in the \"EX-BIOC\" portfolio", () => titledRunSaved(page, 4, "Feed 4", "BDD-ET", "EX-BIOC"));
    await session.step(16, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(17, "And user clicks on scratch evidence project card", () => clickOn(page, el("scratch evidence project card")));
    await session.step(18, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(19, "And user clicks on choose process description run button", () => clickOn(page, el("choose process description run button")));
    await session.step(20, "Then choose run dialog should be visible", () => shouldBe(page, el("choose run dialog"), "visible"));
    await session.step(21, "When user selects \"Feed 4\" in selected run input", () => selectIn(page, "Feed 4", el("selected run input")));
    await session.step(22, "And user enters \"Higher titer at the same temperature\" into selection rationale input", () => enterInto(page, "Higher titer at the same temperature", el("selection rationale input")));
    await session.step(23, "And user clicks on select run button", () => clickOn(page, el("select run button")));
    await session.step(24, "Then choose run dialog should not be visible", () => shouldNotBe(page, el("choose run dialog"), "visible"));
    await session.step(25, "When user clicks on generate process description evidence button", () => clickOn(page, el("generate process description evidence button")));
    await session.step(26, "Then process description document card should contain text \"Evidence table, Locked process description · v1 · draft\"", () => shouldContainText(page, el("process description document card"), "Evidence table, Locked process description · v1 · draft"));
    await session.step(27, "And process description document card should contain text \"Table\"", () => shouldContainText(page, el("process description document card"), "Table"));
  });
  test("Another run is selected, the table goes stale, and is generated again", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(30, "Given a scratch project \"BDD-ET2\" at stage \"D4\" in the \"EX-BIOC\" portfolio", () => scratchProjectAtStage(page, "BDD-ET2", "D4", "EX-BIOC"));
    await session.step(31, "And the workspace of the scratch project \"BDD-ET2\" in the \"EX-BIOC\" portfolio is built", () => scratchWorkspace(page, "BDD-ET2", "EX-BIOC"));
    await session.step(32, "And a run of the example model at feed rate 2 titled \"Feed 2\" saved to the scratch project \"BDD-ET2\" in the \"EX-BIOC\" portfolio", () => titledRunSaved(page, 2, "Feed 2", "BDD-ET2", "EX-BIOC"));
    await session.step(33, "And a run of the example model at feed rate 4 titled \"Feed 4\" saved to the scratch project \"BDD-ET2\" in the \"EX-BIOC\" portfolio", () => titledRunSaved(page, 4, "Feed 4", "BDD-ET2", "EX-BIOC"));
    await session.step(34, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(35, "And user clicks on scratch stale evidence project card", () => clickOn(page, el("scratch stale evidence project card")));
    await session.step(36, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(37, "And user clicks on choose process description run button", () => clickOn(page, el("choose process description run button")));
    await session.step(38, "And user selects \"Feed 4\" in selected run input", () => selectIn(page, "Feed 4", el("selected run input")));
    await session.step(39, "And user enters \"Higher titer at the same temperature\" into selection rationale input", () => enterInto(page, "Higher titer at the same temperature", el("selection rationale input")));
    await session.step(40, "And user clicks on select run button", () => clickOn(page, el("select run button")));
    await session.step(41, "Then choose run dialog should not be visible", () => shouldNotBe(page, el("choose run dialog"), "visible"));
    await session.step(42, "When user clicks on generate process description evidence button", () => clickOn(page, el("generate process description evidence button")));
    await session.step(43, "Then process description document card should contain text \"Evidence table, Locked process description · v1 · draft\"", () => shouldContainText(page, el("process description document card"), "Evidence table, Locked process description · v1 · draft"));
    await session.step(44, "And process description document card should not contain text \"stale\"", () => shouldNotContainText(page, el("process description document card"), "stale"));
    await session.step(45, "When user clicks on choose process description run button", () => clickOn(page, el("choose process description run button")));
    await session.step(46, "And user selects \"Feed 2\" in selected run input", () => selectIn(page, "Feed 2", el("selected run input")));
    await session.step(47, "And user enters \"Lower feed holds the titer\" into selection rationale input", () => enterInto(page, "Lower feed holds the titer", el("selection rationale input")));
    await session.step(48, "And user clicks on select run button", () => clickOn(page, el("select run button")));
    await session.step(49, "Then choose run dialog should not be visible", () => shouldNotBe(page, el("choose run dialog"), "visible"));
    await session.step(50, "And process description document card should contain text \"Evidence changed since written: another run selected: Feed 2 (was Feed 4)\"", () => shouldContainText(page, el("process description document card"), "Evidence changed since written: another run selected: Feed 2 (was Feed 4)"));
    await session.step(51, "When user clicks on regenerate process description evidence link", () => clickOn(page, el("regenerate process description evidence link")));
    await session.step(52, "Then process description document card should contain text \"Evidence table, Locked process description · v2 · draft\"", () => shouldContainText(page, el("process description document card"), "Evidence table, Locked process description · v2 · draft"));
    await session.step(53, "And process description document card should contain text \"Evidence table, Locked process description · v1 · superseded\"", () => shouldContainText(page, el("process description document card"), "Evidence table, Locked process description · v1 · superseded"));
    await session.step(54, "And process description document card should not contain text \"Evidence changed since written\"", () => shouldNotContainText(page, el("process description document card"), "Evidence changed since written"));
  });
});
