/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/find.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addressEndsWith, openApp, recentsCleared} from '../bindings/steps.js';
import {clickOn, enterInto, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Find a program from the landing page", () => {
  const session = feature(test, "features/find.feature", import.meta.url);
  test("Find a program of another portfolio by its target, then open it again from Recent", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given the user's recent programs are cleared", () => recentsCleared(page));
    await session.step(12, "When user enters \"ALP-1\" into program search input", () => enterInto(page, "ALP-1", el("program search input")));
    await session.step(13, "Then search results should contain text \"Antibody Alpha\"", () => shouldContainText(page, el("search results"), "Antibody Alpha"));
    await session.step(14, "And search results should contain text \"target: ALP-1\"", () => shouldContainText(page, el("search results"), "target: ALP-1"));
    await session.step(15, "When user clicks on EXB-201 search hit", () => clickOn(page, el("EXB-201 search hit")));
    await session.step(16, "Then board header should contain text \"EX-BIOC-V1\"", () => shouldContainText(page, el("board header"), "EX-BIOC-V1"));
    await session.step(17, "And program detail should contain text \"Antibody Alpha\"", () => shouldContainText(page, el("program detail"), "Antibody Alpha"));
    await session.step(18, "And the address should end with \"/EX-BIOC/EXB-201\"", () => addressEndsWith(page, "/EX-BIOC/EXB-201"));
    await session.step(19, "When user clicks on portfolios link", () => clickOn(page, el("portfolios link")));
    await session.step(20, "Then recent programs should contain text \"Antibody Alpha\"", () => shouldContainText(page, el("recent programs"), "Antibody Alpha"));
    await session.step(21, "When user clicks on EXB-201 recent", () => clickOn(page, el("EXB-201 recent")));
    await session.step(22, "Then board header should contain text \"EX-BIOC-V1\"", () => shouldContainText(page, el("board header"), "EX-BIOC-V1"));
    await session.step(23, "And program detail should contain text \"Antibody Alpha\"", () => shouldContainText(page, el("program detail"), "Antibody Alpha"));
  });
});
