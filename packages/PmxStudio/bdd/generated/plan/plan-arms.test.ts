/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/plan/plan-arms.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {finishedRunning} from '../../bindings/pmx.js';
import {openApp} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, enterInto, selectIn, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Plan's treatment arms", () => {
  const session = feature(test, "features/plan/plan-arms.feature", import.meta.url);
  test("Plan's treatment arms", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And user opens PMX Studio", () => openApp(page));
    await session.step(9, "Then the study should have finished running", () => finishedRunning(page));
    await run.scenario("One arm needs no correction", async () => {
      await session.step(12, "Then \"Correct for multiple comparisons\" input should be disabled", () => shouldBe(page, el("\"Correct for multiple comparisons\" input"), "disabled"));
      await session.step(13, "And \"Add treatment arm\" button should be enabled", () => shouldBe(page, el("\"Add treatment arm\" button"), "enabled"));
    });
    await run.scenario("Name, therapy, metric and effect reach the plan", async () => {
      await session.step(16, "When user enters \"Compound A\" into \"Name\" input", () => enterInto(page, "Compound A", el("\"Name\" input")));
      await session.step(17, "Then stale note should be visible", () => shouldBe(page, el("stale note"), "visible"));
      await session.step(18, "When user selects \"Cell therapy\" in \"Therapy\" input", () => selectIn(page, "Cell therapy", el("\"Therapy\" input")));
      await session.step(19, "And user clicks on tc choice", () => clickOn(page, el("tc choice")));
      await session.step(20, "And user enters \"40\" into \"Expected\" input", () => enterInto(page, "40", el("\"Expected\" input")));
      await session.step(21, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(22, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(23, "And \"Experimental groups\" section should contain text \"Compound A\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "Compound A"));
      await session.step(24, "And \"Experimental groups\" section should contain text \"T/C 40%\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "T/C 40%"));
      await session.step(25, "When user enters \"99\" into \"Expected\" input", () => enterInto(page, "99", el("\"Expected\" input")));
      await session.step(26, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(27, "When user clears \"Expected\" input", () => clearField(page, el("\"Expected\" input")));
      await session.step(28, "Then form run button should be disabled", () => shouldBe(page, el("form run button"), "disabled"));
      await session.step(29, "When user enters \"40\" into \"Expected\" input", () => enterInto(page, "40", el("\"Expected\" input")));
      await session.step(30, "Then form run button should be enabled", () => shouldBe(page, el("form run button"), "enabled"));
    });
    await run.scenario("Arms are added up to the limit and removed", async () => {
      await session.step(33, "When user clicks on \"Add treatment arm\" button", () => clickOn(page, el("\"Add treatment arm\" button")));
      await session.step(34, "Then \"Correct for multiple comparisons\" input should be enabled", () => shouldBe(page, el("\"Correct for multiple comparisons\" input"), "enabled"));
      await session.step(35, "And form run button should be enabled", () => shouldBe(page, el("form run button"), "enabled"));
      await session.step(36, "And \"Multiple comparisons\" input should contain text \"Dunnett\"", () => shouldContainText(page, el("\"Multiple comparisons\" input"), "Dunnett"));
      await session.step(37, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(38, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(39, "And \"Experimental groups\" section should contain text \"Treatment 2\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "Treatment 2"));
      await session.step(40, "When user clicks on \"Add treatment arm\" button", () => clickOn(page, el("\"Add treatment arm\" button")));
      await session.step(41, "And user clicks on \"Add treatment arm\" button", () => clickOn(page, el("\"Add treatment arm\" button")));
      await session.step(42, "And user clicks on \"Add treatment arm\" button", () => clickOn(page, el("\"Add treatment arm\" button")));
      await session.step(43, "And user clicks on \"Add treatment arm\" button", () => clickOn(page, el("\"Add treatment arm\" button")));
      await session.step(44, "Then \"Add treatment arm\" button should be disabled", () => shouldBe(page, el("\"Add treatment arm\" button"), "disabled"));
      await session.step(45, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(46, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(47, "And \"Experimental groups\" section should contain text \"Treatment 6\"", () => shouldContainText(page, el("\"Experimental groups\" section"), "Treatment 6"));
    });
    await run.scenario("Correction off", async () => {
      await session.step(50, "When user unchecks \"Correct for multiple comparisons\" input", () => uncheck(page, el("\"Correct for multiple comparisons\" input")));
      await session.step(51, "Then \"Multiple comparisons\" input should contain text \"exploratory\"", () => shouldContainText(page, el("\"Multiple comparisons\" input"), "exploratory"));
      await session.step(52, "When user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(53, "Then the study should have finished running", () => finishedRunning(page));
      await session.step(54, "When user clicks on \"Protocol\" tab", () => clickOn(page, el("\"Protocol\" tab")));
      await session.step(55, "Then \"Study protocol\" section should contain text \"None: each of the 6 comparisons\"", () => shouldContainText(page, el("\"Study protocol\" section"), "None: each of the 6 comparisons"));
      await session.step(56, "When user clicks on \"Design\" tab", () => clickOn(page, el("\"Design\" tab")));
      await session.step(57, "And user checks \"Correct for multiple comparisons\" input", () => check(page, el("\"Correct for multiple comparisons\" input")));
      await session.step(58, "And user clicks on form run button", () => clickOn(page, el("form run button")));
      await session.step(59, "Then the study should have finished running", () => finishedRunning(page));
    });
    run.finish();
  });
});
