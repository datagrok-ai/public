/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/provenance.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProjectAtStage} from '../bindings/steps.js';
import {check, clickOn, enterInto, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Provenance", () => {
  const session = feature(test, "features/provenance.feature", import.meta.url);
  test("The provenance of a decision: its readiness and snapshot, the data, and the project it moved", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given a scratch project \"BDD-PROV\" at stage \"TV\" in the \"EX-SMD\" portfolio", () => scratchProjectAtStage(page, "BDD-PROV", "TV", "EX-SMD"));
    await session.step(12, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(13, "And user clicks on provenance project card", () => clickOn(page, el("provenance project card")));
    await session.step(14, "And user clicks on assay criterion met button", () => clickOn(page, el("assay criterion met button")));
    await session.step(15, "And user checks evidence snapshot checkbox", () => check(page, el("evidence snapshot checkbox")));
    await session.step(16, "And user enters \"Assay validated\" into rationale input", () => enterInto(page, "Assay validated", el("rationale input")));
    await session.step(17, "And user clicks on advance button", () => clickOn(page, el("advance button")));
    await session.step(18, "And user clicks on latest decision provenance link", () => clickOn(page, el("latest decision provenance link")));
    await session.step(19, "Then provenance dialog should be visible", () => shouldBe(page, el("provenance dialog"), "visible"));
    await session.step(20, "And provenance start should contain text \"Decision: advance\"", () => shouldContainText(page, el("provenance start"), "Decision: advance"));
    await session.step(21, "And provenance rests on should contain text \"Criterion: Primary assay validated\"", () => shouldContainText(page, el("provenance rests on"), "Criterion: Primary assay validated"));
    await session.step(22, "And provenance rests on should contain text \"Snapshot: BDD-PROV\"", () => shouldContainText(page, el("provenance rests on"), "Snapshot: BDD-PROV"));
    await session.step(23, "And provenance rests on should contain text \"Data: SAR table of BDD-PROV\"", () => shouldContainText(page, el("provenance rests on"), "Data: SAR table of BDD-PROV"));
    await session.step(24, "And provenance used by should contain text \"Program: BDD-PROV\"", () => shouldContainText(page, el("provenance used by"), "Program: BDD-PROV"));
    await session.step(25, "When user clicks on first used-by record", () => clickOn(page, el("first used-by record")));
    await session.step(26, "Then provenance dialog should not be visible", () => shouldNotBe(page, el("provenance dialog"), "visible"));
    await session.step(27, "And program title should contain text \"BDD-PROV\"", () => shouldContainText(page, el("program title"), "BDD-PROV"));
  });
});
