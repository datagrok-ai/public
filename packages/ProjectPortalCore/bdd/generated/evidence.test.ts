/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/evidence.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProjectAtStage, scratchWorkspace} from '../bindings/steps.js';
import {check, clickOn, enterInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Decisions with evidence", () => {
  const session = feature(test, "features/evidence.feature", import.meta.url);
  test("An approver advances a project with its readiness and a link as evidence", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given a scratch project \"BDD-EVIDENCE\" at stage \"TV\" in the \"EX-SMD\" portfolio", () => scratchProjectAtStage(page, "BDD-EVIDENCE", "TV", "EX-SMD"));
    await session.step(12, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(13, "And user clicks on evidence project card", () => clickOn(page, el("evidence project card")));
    await session.step(14, "And user clicks on assay criterion met button", () => clickOn(page, el("assay criterion met button")));
    await session.step(15, "Then decision evidence should contain text \"Readiness of Hit Finding: 1 criterion\"", () => shouldContainText(page, el("decision evidence"), "Readiness of Hit Finding: 1 criterion"));
    await session.step(16, "When user enters \"https://example.org/reviews/hf-gate\" into evidence link input", () => enterInto(page, "https://example.org/reviews/hf-gate", el("evidence link input")));
    await session.step(17, "And user clicks on attach evidence button", () => clickOn(page, el("attach evidence button")));
    await session.step(18, "Then decision evidence should contain text \"Link: example.org/…/hf-gate\"", () => shouldContainText(page, el("decision evidence"), "Link: example.org/…/hf-gate"));
    await session.step(19, "When user enters \"Assay validated, Z' 0.7\" into rationale input", () => enterInto(page, "Assay validated, Z' 0.7", el("rationale input")));
    await session.step(20, "And user clicks on advance button", () => clickOn(page, el("advance button")));
    await session.step(21, "Then latest decision evidence should contain text \"Link: example.org/…/hf-gate\"", () => shouldContainText(page, el("latest decision evidence"), "Link: example.org/…/hf-gate"));
    await session.step(22, "And latest decision evidence should contain text \"Criterion: Primary assay validated\"", () => shouldContainText(page, el("latest decision evidence"), "Criterion: Primary assay validated"));
  });
  test("A decision takes a snapshot of the project now, which opens read-only from the history", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(25, "Given a scratch project \"BDD-EVIDENCE\" at stage \"TV\" in the \"EX-SMD\" portfolio", () => scratchProjectAtStage(page, "BDD-EVIDENCE", "TV", "EX-SMD"));
    await session.step(26, "And the workspace of the scratch project \"BDD-EVIDENCE\" in the \"EX-SMD\" portfolio is built", () => scratchWorkspace(page, "BDD-EVIDENCE", "EX-SMD"));
    await session.step(27, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(28, "And user clicks on evidence project card", () => clickOn(page, el("evidence project card")));
    await session.step(29, "And user clicks on assay criterion met button", () => clickOn(page, el("assay criterion met button")));
    await session.step(30, "And user checks evidence snapshot checkbox", () => check(page, el("evidence snapshot checkbox")));
    await session.step(31, "And user enters \"Assay validated\" into rationale input", () => enterInto(page, "Assay validated", el("rationale input")));
    await session.step(32, "And user clicks on advance button", () => clickOn(page, el("advance button")));
    await session.step(33, "Then latest decision evidence should contain text \"Snapshot: BDD-EVIDENCE\"", () => shouldContainText(page, el("latest decision evidence"), "Snapshot: BDD-EVIDENCE"));
    await session.step(34, "And program snapshots should contain text \"Target Validation\"", () => shouldContainText(page, el("program snapshots"), "Target Validation"));
    await session.step(35, "When user clicks on latest decision snapshot link", () => clickOn(page, el("latest decision snapshot link")));
    await session.step(36, "Then snapshot dialog should be visible", () => shouldBe(page, el("snapshot dialog"), "visible"));
    await session.step(37, "And snapshot stage should contain text \"Target Validation\"", () => shouldContainText(page, el("snapshot stage"), "Target Validation"));
    await session.step(38, "And snapshot criteria should contain text \"Primary assay validated\"", () => shouldContainText(page, el("snapshot criteria"), "Primary assay validated"));
    await session.step(39, "And snapshot table should contain text \"the file holds what was frozen\"", () => shouldContainText(page, el("snapshot table"), "the file holds what was frozen"));
  });
});
