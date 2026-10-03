/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/dmta.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {exampleDesignSets, openApp, openTableStarting} from '../bindings/steps.js';
import {clickOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("DMTA status", () => {
  const session = feature(test, "features/dmta.feature", import.meta.url);
  test("KIN-A's loop and cycle time on EXS-KIN1, and the status per compound", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Given the example design sets are loaded", () => exampleDesignSets(page));
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(15, "Then dmta heading should contain text \"DMTA (18)\"", () => shouldContainText(page, el("dmta heading"), "DMTA (18)"));
    await session.step(16, "And dmta cycles should contain text \"KIN-A: 14 compounds · 9 designed, 6 requested, 8 made, 3 tested\"", () => shouldContainText(page, el("dmta cycles"), "KIN-A: 14 compounds · 9 designed, 6 requested, 8 made, 3 tested"));
    await session.step(17, "And dmta cycles should contain text \"cycle 14 days (median of 3)\"", () => shouldContainText(page, el("dmta cycles"), "cycle 14 days (median of 3)"));
    await session.step(18, "When user clicks on open dmta status button", () => clickOn(page, el("open dmta status button")));
    await session.step(19, "Then the open table should have 18 rows and a name starting with \"EXS-KIN1 DMTA status\"", () => openTableStarting(page, 18, "EXS-KIN1 DMTA status"));
  });
});
