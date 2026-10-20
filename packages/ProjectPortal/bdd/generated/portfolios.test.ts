/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/portfolios.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addressEndsWith, boardStages, currentSpace, openApp} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Several portfolios on one instance", () => {
  const session = feature(test, "features/portfolios.feature", import.meta.url);
  test("Every configured portfolio is listed", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(12, "Then portfolio list should contain text \"Example Small Molecule Discovery\"", () => shouldContainText(page, el("portfolio list"), "Example Small Molecule Discovery"));
    await session.step(13, "And portfolio list should contain text \"Example Biologics Discovery\"", () => shouldContainText(page, el("portfolio list"), "Example Biologics Discovery"));
    await session.step(14, "And portfolio list should contain text \"Example Biologics Development\"", () => shouldContainText(page, el("portfolio list"), "Example Biologics Development"));
  });
  test("Biologics development runs its gates, owned by the registry", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(17, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(18, "Then board header should contain text \"EX-BIOC-V1\"", () => shouldContainText(page, el("board header"), "EX-BIOC-V1"));
    await session.step(19, "And board header should contain text \"stages from registry\"", () => shouldContainText(page, el("board header"), "stages from registry"));
    await session.step(20, "And the board should show the stages \"D1, D2, D3, D4, D5, D6, D7\"", () => boardStages(page, "D1, D2, D3, D4, D5, D6, D7"));
    await session.step(21, "And EXB-201 card in gate D4 column should be visible", () => shouldBe(page, el("EXB-201 card in gate D4 column"), "visible"));
  });
  test("Small molecule discovery runs its own stages, decided in the portal", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(24, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(25, "Then board header should contain text \"EX-SMD-V1\"", () => shouldContainText(page, el("board header"), "EX-SMD-V1"));
    await session.step(26, "And board header should contain text \"stages decided here\"", () => shouldContainText(page, el("board header"), "stages decided here"));
    await session.step(27, "And the board should show the stages \"TV, HF, H2L, LO, CS, PCD\"", () => boardStages(page, "TV, HF, H2L, LO, CS, PCD"));
    await session.step(28, "And EXS-KIN1 card in hit-to-lead column should be visible", () => shouldBe(page, el("EXS-KIN1 card in hit-to-lead column"), "visible"));
  });
  test("Biologics discovery reads the shared registry but not its gates", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(31, "When user clicks on biologics discovery portfolio", () => clickOn(page, el("biologics discovery portfolio")));
    await session.step(32, "Then board header should contain text \"EX-BIOD-V1\"", () => shouldContainText(page, el("board header"), "EX-BIOD-V1"));
    await session.step(33, "And the board should show the stages \"TI, LG, LO, DA, CN\"", () => boardStages(page, "TI, LG, LO, DA, CN"));
    await session.step(34, "And EXB-101 card in lead generation column should be visible", () => shouldBe(page, el("EXB-101 card in lead generation column"), "visible"));
    await session.step(35, "When user clicks on EXB-101 card", () => clickOn(page, el("EXB-101 card")));
    await session.step(36, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(37, "Then program detail should contain text \"Registry stage: 3\"", () => shouldContainText(page, el("program detail"), "Registry stage: 3"));
  });
  test("A portfolio has its own address, and a program opens its workspace", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(40, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(41, "Then the address should end with \"/Portal/EX-SMD\"", () => addressEndsWith(page, "/Portal/EX-SMD"));
    await session.step(42, "When user clicks on portfolios link", () => clickOn(page, el("portfolios link")));
    await session.step(43, "And user clicks on biologics discovery portfolio", () => clickOn(page, el("biologics discovery portfolio")));
    await session.step(44, "Then the address should end with \"/Portal/EX-BIOD\"", () => addressEndsWith(page, "/Portal/EX-BIOD"));
    await session.step(45, "When user clicks on portfolios link", () => clickOn(page, el("portfolios link")));
    await session.step(46, "And user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(47, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(48, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(49, "And user clicks on open workspace button", () => clickOn(page, el("open workspace button")));
    await session.step(50, "Then the open view should be the space \"EXS-KIN1\"", () => currentSpace(page, "EXS-KIN1"));
  });
  test("A program shows the work units its portfolio uses", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(53, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(54, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(55, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(56, "Then work units should contain text \"Campaigns (2)\"", () => shouldContainText(page, el("work units"), "Campaigns (2)"));
    await session.step(57, "And KIN-A series should contain text \"validated\"", () => shouldContainText(page, el("KIN-A series"), "validated"));
    await session.step(58, "When user clicks on portfolios link", () => clickOn(page, el("portfolios link")));
    await session.step(59, "And user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(60, "And user clicks on EXB-201 card", () => clickOn(page, el("EXB-201 card")));
    await session.step(61, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(62, "Then work units should contain text \"Process stages (4)\"", () => shouldContainText(page, el("work units"), "Process stages (4)"));
    await session.step(63, "And TOX-201 study should contain text \"completed\"", () => shouldContainText(page, el("TOX-201 study"), "completed"));
  });
  test("A project lists its compounds by role, never as folders", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(66, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(67, "And user clicks on EXS-KIN1 card", () => clickOn(page, el("EXS-KIN1 card")));
    await session.step(68, "And user clicks on open program button", () => clickOn(page, el("open program button")));
    await session.step(69, "Then compound roles should contain text \"6 lead\"", () => shouldContainText(page, el("compound roles"), "6 lead"));
    await session.step(70, "And EXC-3001 compound should contain text \"tool\"", () => shouldContainText(page, el("EXC-3001 compound"), "tool"));
  });
  test("No data access notice when every group can query the portfolio's data", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(73, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(74, "Then board header should contain text \"EX-SMD-V1\"", () => shouldContainText(page, el("board header"), "EX-SMD-V1"));
    await session.step(75, "And data access notice should not be visible", () => shouldNotBe(page, el("data access notice"), "visible"));
  });
});
