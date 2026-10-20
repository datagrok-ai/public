/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/predictions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {notScored, openApp, openTableEnding} from '../bindings/steps.js';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Predictions", () => {
  const session = feature(test, "features/predictions.feature", import.meta.url);
  test("Score a project's compounds and set the predictions against the measurements", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given the example solubility model has not scored \"EXS-KIN1\" in the \"EX-SMD\" portfolio", () => notScored(page, "EXS-KIN1", "EX-SMD"));
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(15, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(16, "And user clicks on LO step", () => clickOn(page, el("LO step")));
    await session.step(17, "And user clicks on score solubility estimate button", () => clickOn(page, el("score solubility estimate button")));
    await session.step(18, "Then score result should contain text \"stored 11 solubility_um\"", () => shouldContainText(page, el("score result"), "stored 11 solubility_um"));
    await session.step(19, "And solubility predictions fit should contain text \"7 of 11 measured\"", () => shouldContainText(page, el("solubility predictions fit"), "7 of 11 measured"));
    await session.step(20, "And missing solubility prediction should contain text \"EXC-1006 34.1 uM\"", () => shouldContainText(page, el("missing solubility prediction"), "EXC-1006 34.1 uM"));
    await session.step(21, "When user clicks on solubility predicted vs measured button", () => clickOn(page, el("solubility predicted vs measured button")));
    await session.step(22, "Then the open table should have 11 rows and a name ending with \"solubility_um: predicted vs measured\"", () => openTableEnding(page, 11, "solubility_um: predicted vs measured"));
  });
});
