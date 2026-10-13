/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/documents.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchFolderFile, scratchProjectAtStage, scratchWorkspace} from '../bindings/steps.js';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Documents over files", () => {
  const session = feature(test, "features/documents.feature", import.meta.url);
  test("A file of the project folder is linked to the process description and approved", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given a scratch project \"BDD-DOCS\" at stage \"D4\" in the \"EX-BIOC\" portfolio", () => scratchProjectAtStage(page, "BDD-DOCS", "D4", "EX-BIOC"));
    await session.step(13, "And the workspace of the scratch project \"BDD-DOCS\" in the \"EX-BIOC\" portfolio is built", () => scratchWorkspace(page, "BDD-DOCS", "EX-BIOC"));
    await session.step(14, "And a file \"process-description-v1.md\" in the folder of the scratch project \"BDD-DOCS\" in the \"EX-BIOC\" portfolio", () => scratchFolderFile(page, "process-description-v1.md", "BDD-DOCS", "EX-BIOC"));
    await session.step(15, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(16, "And user clicks on scratch docs project card", () => clickOn(page, el("scratch docs project card")));
    await session.step(17, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(18, "Then process description document card should contain text \"no document yet\"", () => shouldContainText(page, el("process description document card"), "no document yet"));
    await session.step(19, "When user clicks on link process description document button", () => clickOn(page, el("link process description document button")));
    await session.step(20, "Then link document dialog should be visible", () => shouldBe(page, el("link document dialog"), "visible"));
    await session.step(21, "When user selects \"process-description-v1.md\" in document file input", () => selectIn(page, "process-description-v1.md", el("document file input")));
    await session.step(22, "And user enters \"Process description, BDD\" into document title input", () => enterInto(page, "Process description, BDD", el("document title input")));
    await session.step(23, "And user enters \"1.0\" into document version input", () => enterInto(page, "1.0", el("document version input")));
    await session.step(24, "And user clicks on link document button", () => clickOn(page, el("link document button")));
    await session.step(25, "Then link document dialog should not be visible", () => shouldNotBe(page, el("link document dialog"), "visible"));
    await session.step(26, "And process description document card should contain text \"Process description, BDD · v1.0 · draft\"", () => shouldContainText(page, el("process description document card"), "Process description, BDD · v1.0 · draft"));
    await session.step(27, "And process description document card should contain text \"in progress\"", () => shouldContainText(page, el("process description document card"), "in progress"));
    await session.step(28, "When user clicks on approve process description document link", () => clickOn(page, el("approve process description document link")));
    await session.step(29, "Then process description document card should contain text \"Process description, BDD · v1.0 · approved\"", () => shouldContainText(page, el("process description document card"), "Process description, BDD · v1.0 · approved"));
    await session.step(30, "And process description document card should contain text \"met\"", () => shouldContainText(page, el("process description document card"), "met"));
  });
});
