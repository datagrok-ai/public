/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/design-sets.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {exampleDesignSets, openApp} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Design sets and selection sets", () => {
  const session = feature(test, "features/design-sets.feature", import.meta.url);
  test("A series' design set opens with its proposals ranked and the chosen ones ticked", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given the example design sets are loaded", () => exampleDesignSets(page));
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(15, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(16, "Then design sets should contain text \"KIN-A para-substituents\"", () => shouldContainText(page, el("design sets"), "KIN-A para-substituents"));
    await session.step(17, "And design sets should contain text \"6 proposed · 4 chosen for synthesis · best KIN-A-D01\"", () => shouldContainText(page, el("design sets"), "6 proposed · 4 chosen for synthesis · best KIN-A-D01"));
    await session.step(18, "And selection sets should contain text \"K1 screen: indazole hits\"", () => shouldContainText(page, el("selection sets"), "K1 screen: indazole hits"));
    await session.step(19, "When user clicks on para design set link", () => clickOn(page, el("para design set link")));
    await session.step(20, "Then design set dialog should be visible", () => shouldBe(page, el("design set dialog"), "visible"));
    await session.step(21, "And design proposals should contain text \"KIN-A-D06\"", () => shouldContainText(page, el("design proposals"), "KIN-A-D06"));
    await session.step(22, "And design set summary should contain text \"scored by predicted pIC50\"", () => shouldContainText(page, el("design set summary"), "scored by predicted pIC50"));
  });
});
