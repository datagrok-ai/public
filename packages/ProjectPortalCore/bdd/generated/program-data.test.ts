/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/program-data.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, openTableEnding, openTableStarting, openTableViewers} from '../bindings/steps.js';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A program's governed data one click away", () => {
  const session = feature(test, "features/program-data.feature", import.meta.url);
  test("The SAR table of a discovery project opens as a table named after it", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(13, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(14, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(15, "Then open data button should be visible", () => shouldBe(page, el("open data button"), "visible"));
    await session.step(16, "When user clicks on open data button", () => clickOn(page, el("open data button")));
    await session.step(17, "Then the open table should have 11 rows and a name starting with \"EXS-KIN1\"", () => openTableStarting(page, 11, "EXS-KIN1"));
    await session.step(18, "And the open table should show the viewer types \"Scatter plot, Histogram\"", () => openTableViewers(page, "Scatter plot, Histogram"));
  });
  test("The assay results open beside it", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(21, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(22, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(23, "And user clicks on open assay results button", () => clickOn(page, el("open assay results button")));
    await session.step(24, "Then the open table should have 8 rows and a name ending with \"Assay results\"", () => openTableEnding(page, 8, "Assay results"));
    await session.step(25, "And the open table should show the viewer types \"Box plot\"", () => openTableViewers(page, "Box plot"));
  });
});
