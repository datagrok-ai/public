/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/templates/triage-forms.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {descriptorsUnlisted, editorShowsError, editorShowsNoError, templatesOnServer} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, selectIn, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openApp, packageInstalled, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hit Triage and PepTriage forms: START and CREATE refuse what they cannot use", () => {
  const session = feature(test, "features/templates/triage-forms.feature", import.meta.url);
  test("Hit Triage and PepTriage forms: START and CREATE refuse what they cannot use", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the descriptors cannot be listed, put back at feature end", () => descriptorsUnlisted(page));
    await session.step(18, "And the \"Hit Triage\" templates are on the server:", () => templatesOnServer(page, "Hit Triage", [["name","key","fields"],[session.text("BDD HT Form {time}"),session.text("BDDHTF{time}"),"Owner: String required"]]), [["name","key","fields"],[session.text("BDD HT Form {time}"),session.text("BDDHTF{time}"),"Owner: String required"]]);
    await session.step(21, "And user opens the \"Hit Triage\" app", () => openApp(page, "Hit Triage"));
    await run.scenario("While a required field is empty, START is disabled", async () => {
      await session.step(24, "When user selects \"BDD HT Form {time}\" in Template input", () => selectIn(page, session.text("BDD HT Form {time}"), el("Template input")));
      await session.step(25, "Then \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(26, "And Start button should be disabled", () => shouldBe(page, el("Start button"), "disabled"));
    });
    await run.scenario("The editor shows its error while the name is empty or an existing template's", async () => {
      await session.step(29, "When user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(30, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(31, "When user types \"BDD HT Form {time}\" into template name input", () => typeInto(page, session.text("BDD HT Form {time}"), el("template name input")));
      await session.step(32, "Then the template editor should show the error \"Template name is empty or already exists\"", () => editorShowsError(page, "Template name is empty or already exists"));
      await session.step(33, "When user types \"BDD HT Free {time}\" into template name input", () => typeInto(page, session.text("BDD HT Free {time}"), el("template name input")));
      await session.step(34, "Then the template editor should not show the error \"Template name is empty or already exists\"", () => editorShowsNoError(page, "Template name is empty or already exists"));
      await session.step(35, "When user clears template name input", () => clearField(page, el("template name input")));
      await session.step(36, "Then the template editor should show the error \"Template name is empty or already exists\"", () => editorShowsError(page, "Template name is empty or already exists"));
      await session.step(37, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(38, "Then the \"Hit Triage\" view should be current", () => viewIsCurrent(page, "Hit Triage"));
    });
    await run.scenario("CREATE is disabled while the name is empty or an existing template's", async () => {
      await session.step(41, "When user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(42, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(43, "And template name input should have the value \"\"", () => shouldHaveValue(page, el("template name input"), ""));
      await session.step(44, "And Create button should be disabled", () => shouldBe(page, el("Create button"), "disabled"));
      await session.step(45, "When user types \"BDD HT Form {time}\" into template name input", () => typeInto(page, session.text("BDD HT Form {time}"), el("template name input")));
      await session.step(46, "Then the template editor should show the error \"Template name is empty or already exists\"", () => editorShowsError(page, "Template name is empty or already exists"));
      await session.step(47, "And Create button should be disabled", () => shouldBe(page, el("Create button"), "disabled"));
      await session.step(48, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(49, "Then the \"Hit Triage\" view should be current", () => viewIsCurrent(page, "Hit Triage"));
    });
    await run.scenario("In PepTriage, CREATE is disabled while the name is empty", async () => {
      await session.step(52, "Given the \"Bio\" package is installed", () => packageInstalled(page, "Bio"));
      await session.step(53, "When user opens the \"PepTriage\" app", () => openApp(page, "PepTriage"));
      await session.step(54, "And user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(55, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(56, "And template name input should have the value \"\"", () => shouldHaveValue(page, el("template name input"), ""));
      await session.step(57, "And Create button should be disabled", () => shouldBe(page, el("Create button"), "disabled"));
      await session.step(58, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(59, "Then the \"PepTriage\" view should be current", () => viewIsCurrent(page, "PepTriage"));
    });
    run.finish();
  });
});
