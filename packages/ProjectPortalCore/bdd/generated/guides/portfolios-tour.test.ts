/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/portfolios-tour.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, scratchProject, sidePanelsHidden} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Several portfolios on one instance", () => {
  const session = feature(test, "features/guides/portfolios-tour.feature", import.meta.url);
  test("Three portfolios, three stage models, two kinds of gate authority", {tag: ["@guide"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And a scratch project \"KINASE-X\" in the \"EX-SMD\" portfolio", () => scratchProject(page, "KINASE-X", "EX-SMD"));
    await session.step(14, "And user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(15, "And the side panels are hidden", () => sidePanelsHidden(page));
    await session.step(16, "Then portfolio list should contain text \"Example Biologics Development\"", () => shouldContainText(page, el("portfolio list"), "Example Biologics Development"));
    await session.step(17, "When user clicks on biologics development portfolio", () => clickOn(page, el("biologics development portfolio")));
    await session.step(18, "Then board header should contain text \"stages from registry\"", () => shouldContainText(page, el("board header"), "stages from registry"));
    await session.step(19, "When user clicks on EXB-201 card", () => clickOn(page, el("EXB-201 card")));
    await session.step(20, "Then gate reasons should contain text \"The registry owns this portfolio's stages\"", () => shouldContainText(page, el("gate reasons"), "The registry owns this portfolio's stages"));
    await session.step(21, "When user clicks on biologics discovery portfolio", () => clickOn(page, el("biologics discovery portfolio")));
    await session.step(22, "Then board header should contain text \"EX-BIOD-V1\"", () => shouldContainText(page, el("board header"), "EX-BIOD-V1"));
    await session.step(23, "When user clicks on EXB-101 card", () => clickOn(page, el("EXB-101 card")));
    await session.step(24, "Then program detail should contain text \"Registry stage: 3\"", () => shouldContainText(page, el("program detail"), "Registry stage: 3"));
    await session.step(25, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(26, "Then board header should contain text \"stages decided here\"", () => shouldContainText(page, el("board header"), "stages decided here"));
    await session.step(27, "When user clicks on KINASE-X card", () => clickOn(page, el("KINASE-X card")));
    await session.step(28, "And user clicks on advance button", () => clickOn(page, el("advance button")));
    await session.step(29, "Then KINASE-X card in target validation column should be visible", () => shouldBe(page, el("KINASE-X card in target validation column"), "visible"));
    await session.step(30, "And gate reasons should contain text \"required criteria of Hit Finding\"", () => shouldContainText(page, el("gate reasons"), "required criteria of Hit Finding"));
    await session.step(31, "When user clicks on assay criterion met button", () => clickOn(page, el("assay criterion met button")));
    await session.step(32, "And user clicks on advance button", () => clickOn(page, el("advance button")));
    await session.step(33, "Then transition history should contain text \"Hit Finding\"", () => shouldContainText(page, el("transition history"), "Hit Finding"));
  });
});
