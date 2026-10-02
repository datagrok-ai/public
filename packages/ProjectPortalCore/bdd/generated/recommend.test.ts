/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/recommend.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp} from '../bindings/steps.js';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Suggested models and missing assays", () => {
  const session = feature(test, "features/recommend.feature", import.meta.url);
  test("What Lead Optimization still needs, and a model that stands in for it", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(13, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(14, "And user clicks on LO step", () => clickOn(page, el("LO step")));
    await session.step(15, "Then missing assays heading should contain text \"Missing assays for Lead Optimization (3)\"", () => shouldContainText(page, el("missing assays heading"), "Missing assays for Lead Optimization (3)"));
    await session.step(16, "And missing solubility should contain text \"solubility_um on 1 row (EXC-1006)\"", () => shouldContainText(page, el("missing solubility"), "solubility_um on 1 row (EXC-1006)"));
    await session.step(17, "And missing solubility should contain text \"Solubility estimate (example)\"", () => shouldContainText(page, el("missing solubility"), "Solubility estimate (example)"));
    await session.step(18, "And suggested models heading should contain text \"Suggested models (2)\"", () => shouldContainText(page, el("suggested models heading"), "Suggested models (2)"));
    await session.step(19, "And suggested solubility estimate should contain text \"predicts solubility_um, missing for the gate\"", () => shouldContainText(page, el("suggested solubility estimate"), "predicts solubility_um, missing for the gate"));
    await session.step(20, "When user clicks on solubility estimate card link", () => clickOn(page, el("solubility estimate card link")));
    await session.step(21, "Then model card question should contain text \"Roughly what solubility\"", () => shouldContainText(page, el("model card question"), "Roughly what solubility"));
  });
});
