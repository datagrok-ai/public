/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/templates/template-picker.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {confirmDialog, descriptorsUnlisted, editorFields, editorStages, templateHas, templateNotOnServer, templatesOnServer, templatesRemoved} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotHaveValue, shouldNotOfferChoice, shouldOfferChoice, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, openApp, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloon, noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hit Design template picker: choosing, adding, cloning and deleting templates", () => {
  const session = feature(test, "features/templates/template-picker.feature", import.meta.url);
  test("Hit Design template picker: choosing, adding, cloning and deleting templates", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the descriptors cannot be listed, put back at feature end", () => descriptorsUnlisted(page));
    await session.step(20, "And the \"Hit Design\" templates are on the server:", () => templatesOnServer(page, "Hit Design", [["name","key","fields","stages","submit","layout"],[session.text("BDD Picker Source {time}"),session.text("BDDPKS{time}"),"Owner: String required","Design, Make","Demo File Submit","bar chart of Stage"],[session.text("BDD Picker Doomed {time}"),session.text("BDDPKD{time}"),"","Design","",""]]), [["name","key","fields","stages","submit","layout"],[session.text("BDD Picker Source {time}"),session.text("BDDPKS{time}"),"Owner: String required","Design, Make","Demo File Submit","bar chart of Stage"],[session.text("BDD Picker Doomed {time}"),session.text("BDDPKD{time}"),"","Design","",""]]);
    await session.step(24, "And the \"Hit Design\" templates named \"BDD Picker Clone {time}\" are removed now and at feature end", () => templatesRemoved(page, "Hit Design", session.text("BDD Picker Clone {time}")));
    await session.step(25, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
    await run.scenario("The picker lists the templates, offers New, Clone and Delete, and a choice builds its form", async () => {
      await session.step(28, "Then Template input should offer the choice \"BDD Picker Source {time}\"", () => shouldOfferChoice(page, el("Template input"), session.text("BDD Picker Source {time}")));
      await session.step(29, "And Template input should offer the choice \"BDD Picker Doomed {time}\"", () => shouldOfferChoice(page, el("Template input"), session.text("BDD Picker Doomed {time}")));
      await session.step(30, "And \"New Template\" icon should be visible", () => shouldBe(page, el("\"New Template\" icon"), "visible"));
      await session.step(31, "And \"Clone template\" icon should be visible", () => shouldBe(page, el("\"Clone template\" icon"), "visible"));
      await session.step(32, "And \"Delete template\" icon should be visible", () => shouldBe(page, el("\"Delete template\" icon"), "visible"));
      await session.step(33, "When user selects \"BDD Picker Source {time}\" in Template input", () => selectIn(page, session.text("BDD Picker Source {time}"), el("Template input")));
      await session.step(34, "Then \"Campaign Name\" input should be visible", () => shouldBe(page, el("\"Campaign Name\" input"), "visible"));
      await session.step(35, "And \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(36, "And Start button should be visible", () => shouldBe(page, el("Start button"), "visible"));
    });
    await run.scenario("The \"+\" icon opens an empty New Template editor, and CANCEL comes back", async () => {
      await session.step(39, "When user closes all views", () => closeAllViews(page));
      await session.step(40, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(41, "And user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(42, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(43, "And the page address should contain \"/new-template\"", () => urlShouldContain(page, "/new-template"));
      await session.step(44, "And template name input should have the value \"\"", () => shouldHaveValue(page, el("template name input"), ""));
      await session.step(45, "And template key input should have the value \"\"", () => shouldHaveValue(page, el("template key input"), ""));
      await session.step(46, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(47, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
    });
    await run.scenario("Clone fills the editor with the selected template's fields and stages", async () => {
      await session.step(50, "When user closes all views", () => closeAllViews(page));
      await session.step(51, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(52, "And user selects \"BDD Picker Source {time}\" in Template input", () => selectIn(page, session.text("BDD Picker Source {time}"), el("Template input")));
      await session.step(53, "Then \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(54, "When user clicks on \"Clone template\" icon", () => clickOn(page, el("\"Clone template\" icon")));
      await session.step(55, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(56, "And the template editor should list the fields \"Owner (String, required)\"", () => editorFields(page, "Owner (String, required)"));
      await session.step(57, "And the template editor should list the stages \"Design, Make\"", () => editorStages(page, "Design, Make"));
      await session.step(58, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(59, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
    });
    await run.scenario("A clone proposes a name and a key no template has", async () => {
      await session.step(62, "When user closes all views", () => closeAllViews(page));
      await session.step(63, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(64, "And user selects \"BDD Picker Source {time}\" in Template input", () => selectIn(page, session.text("BDD Picker Source {time}"), el("Template input")));
      await session.step(65, "Then \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(66, "When user clicks on \"Clone template\" icon", () => clickOn(page, el("\"Clone template\" icon")));
      await session.step(67, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(68, "And template name input should have the value \"BDD Picker Source {time} (copy)\"", () => shouldHaveValue(page, el("template name input"), session.text("BDD Picker Source {time} (copy)")));
      await session.step(69, "And template key input should not have the value \"BDDPKS{time}\"", () => shouldNotHaveValue(page, el("template key input"), session.text("BDDPKS{time}")));
      await session.step(70, "And Create button should be enabled", () => shouldBe(page, el("Create button"), "enabled"));
      await session.step(71, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(72, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
    });
    await run.scenario("A clone keeps the template's submit function", async () => {
      await session.step(75, "When user closes all views", () => closeAllViews(page));
      await session.step(76, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(77, "And user selects \"BDD Picker Source {time}\" in Template input", () => selectIn(page, session.text("BDD Picker Source {time}"), el("Template input")));
      await session.step(78, "Then \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(79, "When user clicks on \"Clone template\" icon", () => clickOn(page, el("\"Clone template\" icon")));
      await session.step(80, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(81, "And \"Submit function\" input should have the value \"Demo File Submit\"", () => shouldHaveValue(page, el("\"Submit function\" input"), "Demo File Submit"));
      await session.step(82, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(83, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
    });
    await run.scenario("A clone saved under a new name carries the template's layout, and CREATE comes back with it selected", async () => {
      await session.step(86, "When user closes all views", () => closeAllViews(page));
      await session.step(87, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(88, "And user selects \"BDD Picker Source {time}\" in Template input", () => selectIn(page, session.text("BDD Picker Source {time}"), el("Template input")));
      await session.step(89, "Then \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(90, "When user clicks on \"Clone template\" icon", () => clickOn(page, el("\"Clone template\" icon")));
      await session.step(91, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(92, "When user types \"BDD Picker Clone {time}\" into template name input", () => typeInto(page, session.text("BDD Picker Clone {time}"), el("template name input")));
      await session.step(93, "And user types \"BDDPKC{time}\" into template key input", () => typeInto(page, session.text("BDDPKC{time}"), el("template key input")));
      await session.step(94, "And user clicks on Create button", () => clickOn(page, el("Create button")));
      await session.step(95, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(96, "And Template input should have the value \"BDD Picker Clone {time}\"", () => shouldHaveValue(page, el("Template input"), session.text("BDD Picker Clone {time}")));
      await session.step(97, "And the \"Hit Design\" template \"BDD Picker Clone {time}\" on the server should have the fields \"Owner (String, required)\"", () => templateHas(page, "Hit Design", session.text("BDD Picker Clone {time}"), "fields", "Owner (String, required)"));
      await session.step(98, "And the \"Hit Design\" template \"BDD Picker Clone {time}\" on the server should have the layout \"a layout\"", () => templateHas(page, "Hit Design", session.text("BDD Picker Clone {time}"), "layout", "a layout"));
    });
    await run.scenario("Delete asks for a confirmation, removes the template and moves the picker to another one, with no error", async () => {
      await session.step(101, "When user closes all views", () => closeAllViews(page));
      await session.step(102, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(103, "And user selects \"BDD Picker Doomed {time}\" in Template input", () => selectIn(page, session.text("BDD Picker Doomed {time}"), el("Template input")));
      await session.step(104, "And user clicks on \"Delete template\" icon", () => clickOn(page, el("\"Delete template\" icon")));
      await session.step(105, "Then \"Delete template\" dialog should contain text \"Are you sure you want to delete template\"", () => shouldContainText(page, el("\"Delete template\" dialog"), "Are you sure you want to delete template"));
      await session.step(106, "And \"Delete template\" dialog should contain text \"BDD Picker Doomed {time}\"", () => shouldContainText(page, el("\"Delete template\" dialog"), session.text("BDD Picker Doomed {time}")));
      await session.step(107, "When user confirms the \"Delete template\" dialog", () => confirmDialog(page, el("the \"Delete template\" dialog")));
      await session.step(108, "Then the \"Hit Design\" template \"BDD Picker Doomed {time}\" should not be on the server", () => templateNotOnServer(page, "Hit Design", session.text("BDD Picker Doomed {time}")));
      await session.step(109, "And an info balloon should have been shown", () => infoBalloon(page));
      await session.step(110, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(111, "And Template input should not offer the choice \"BDD Picker Doomed {time}\"", () => shouldNotOfferChoice(page, el("Template input"), session.text("BDD Picker Doomed {time}")));
      await session.step(112, "And Template input should not have the value \"\"", () => shouldNotHaveValue(page, el("Template input"), ""));
    });
    run.finish();
  });
});
