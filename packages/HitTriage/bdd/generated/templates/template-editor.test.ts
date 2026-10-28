/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/templates/template-editor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addRequiredField, addStage, descriptorsUnlisted, editorFields, editorShowsError, editorShowsNoError, editorStages, templateHas, templateNotOnServer, templatesOnServer, templatesRemoved, typePendingField, typePendingStage, viewNotOpen} from '../../bindings/landing.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, selectIn, shouldBe, shouldHaveValue, shouldNotOfferChoice, shouldOfferChoice, typeInto, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {closeAllViews, currentViewShows, openApp, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Hit Design template editor: name and key checks, fields, stages, submit function, layout, CREATE and CANCEL", () => {
  const session = feature(test, "features/templates/template-editor.feature", import.meta.url);
  test("Hit Design template editor: name and key checks, fields, stages, submit function, layout, CREATE and CANCEL", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the descriptors cannot be listed, put back at feature end", () => descriptorsUnlisted(page));
    await session.step(21, "And the \"Hit Design\" templates are on the server:", () => templatesOnServer(page, "Hit Design", [["name","key"],[session.text("BDD Editor Existing {time}"),session.text("BDDEDX{time}")]]), [["name","key"],[session.text("BDD Editor Existing {time}"),session.text("BDDEDX{time}")]]);
    await session.step(24, "And the \"Hit Design\" templates named \"BDD Editor New {time}, BDD Editor Cancelled {time}\" are removed now and at feature end", () => templatesRemoved(page, "Hit Design", session.text("BDD Editor New {time}, BDD Editor Cancelled {time}")));
    await session.step(25, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
    await run.scenario("The editor has Details, Additional fields, Stages, Compute and Submit, with CANCEL and CREATE", async () => {
      await session.step(28, "When user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(29, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(30, "And \"Hit Design\" link should be visible", () => shouldBe(page, el("\"Hit Design\" link"), "visible"));
      await session.step(31, "And template name input should be visible", () => shouldBe(page, el("template name input"), "visible"));
      await session.step(32, "And template key input should be visible", () => shouldBe(page, el("template key input"), "visible"));
      await session.step(33, "And Layout input should be visible", () => shouldBe(page, el("Layout input"), "visible"));
      await session.step(34, "And \"Submit function\" input should be visible", () => shouldBe(page, el("\"Submit function\" input"), "visible"));
      await session.step(35, "And the current view should show the text \"Additional fields\"", () => currentViewShows(page, "Additional fields"));
      await session.step(36, "And the current view should show the text \"Stages\"", () => currentViewShows(page, "Stages"));
      await session.step(37, "And the current view should show the text \"Compute\"", () => currentViewShows(page, "Compute"));
      await session.step(38, "And the current view should show the text \"Could not list the descriptors: the descriptors are not available\"", () => currentViewShows(page, "Could not list the descriptors: the descriptors are not available"));
      await session.step(39, "And Cancel button should be visible", () => shouldBe(page, el("Cancel button"), "visible"));
      await session.step(40, "And Create button should be visible", () => shouldBe(page, el("Create button"), "visible"));
      await session.step(41, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(42, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
    });
    await run.scenario("A name or a key that is empty or taken shows its error, a free one clears it", async () => {
      await session.step(45, "When user closes all views", () => closeAllViews(page));
      await session.step(46, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(47, "And user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(48, "And user types \"BDD Editor Existing {time}\" into template name input", () => typeInto(page, session.text("BDD Editor Existing {time}"), el("template name input")));
      await session.step(49, "Then the template editor should show the error \"A template with this name already exists\"", () => editorShowsError(page, "A template with this name already exists"));
      await session.step(50, "When user types \"BDD Editor Free {time}\" into template name input", () => typeInto(page, session.text("BDD Editor Free {time}"), el("template name input")));
      await session.step(51, "Then the template editor should not show the error \"A template with this name already exists\"", () => editorShowsNoError(page, "A template with this name already exists"));
      await session.step(52, "When user types \"BDDEDX{time}\" into template key input", () => typeInto(page, session.text("BDDEDX{time}"), el("template key input")));
      await session.step(53, "Then the template editor should show the error \"BDDEDX{time} is used by BDD Editor Existing {time}\"", () => editorShowsError(page, session.text("BDDEDX{time} is used by BDD Editor Existing {time}")));
      await session.step(54, "When user types \"BDDEDF{time}\" into template key input", () => typeInto(page, session.text("BDDEDF{time}"), el("template key input")));
      await session.step(55, "Then the template editor should not show the error \"BDDEDX{time} is used by BDD Editor Existing {time}\"", () => editorShowsNoError(page, session.text("BDDEDX{time} is used by BDD Editor Existing {time}")));
      await session.step(56, "When user clears template name input", () => clearField(page, el("template name input")));
      await session.step(57, "Then the template editor should show the error \"Value can't be empty\"", () => editorShowsError(page, "Value can't be empty"));
      await session.step(58, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(59, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
    });
    await run.scenario("While the name is taken, CREATE is disabled", async () => {
      await session.step(62, "When user closes all views", () => closeAllViews(page));
      await session.step(63, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(64, "And user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(65, "Then Create button should be disabled", () => shouldBe(page, el("Create button"), "disabled"));
      await session.step(66, "When user types \"BDD Editor Existing {time}\" into template name input", () => typeInto(page, session.text("BDD Editor Existing {time}"), el("template name input")));
      await session.step(67, "And user types \"BDDEDF{time}\" into template key input", () => typeInto(page, session.text("BDDEDF{time}"), el("template key input")));
      await session.step(68, "Then the template editor should show the error \"A template with this name already exists\"", () => editorShowsError(page, "A template with this name already exists"));
      await session.step(69, "And Create button should be disabled", () => shouldBe(page, el("Create button"), "disabled"));
      await session.step(70, "When user types \"BDD Editor Free {time}\" into template name input", () => typeInto(page, session.text("BDD Editor Free {time}"), el("template name input")));
      await session.step(71, "Then Create button should be enabled", () => shouldBe(page, el("Create button"), "enabled"));
      await session.step(72, "When user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(73, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
    });
    await run.scenario("CREATE writes the fields, stages, submit function and layout, the rows being typed included, and comes back with the template selected", async () => {
      await session.step(76, "When user closes all views", () => closeAllViews(page));
      await session.step(77, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(78, "And user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(79, "And user types \"BDD Editor New {time}\" into template name input", () => typeInto(page, session.text("BDD Editor New {time}"), el("template name input")));
      await session.step(80, "And user types \"BDDEDN{time}\" into template key input", () => typeInto(page, session.text("BDDEDN{time}"), el("template key input")));
      await session.step(81, "And user uploads \"fixtures/bar-chart.layout\" through \"Open file\" icon in Layout input", () => uploadThrough(page, "fixtures/bar-chart.layout", el("\"Open file\" icon in Layout input")));
      await session.step(82, "And user adds the required field \"Owner\" of type \"String\" in the template editor", () => addRequiredField(page, "Owner", "String"));
      await session.step(83, "Then the template editor should list the fields \"Owner (String, required)\"", () => editorFields(page, "Owner (String, required)"));
      await session.step(84, "When user types the field \"Due\" of type \"Date\" in the template editor without adding it", () => typePendingField(page, "Due", "Date"));
      await session.step(85, "And user adds the stage \"Design\" in the template editor", () => addStage(page, "Design"));
      await session.step(86, "Then the template editor should list the stages \"Design\"", () => editorStages(page, "Design"));
      await session.step(87, "When user types the stage \"Make\" in the template editor without adding it", () => typePendingStage(page, "Make"));
      await session.step(88, "And user selects \"Demo File Submit\" in \"Submit function\" input", () => selectIn(page, "Demo File Submit", el("\"Submit function\" input")));
      await session.step(89, "And user clicks on Create button", () => clickOn(page, el("Create button")));
      await session.step(90, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(91, "And Template input should have the value \"BDD Editor New {time}\"", () => shouldHaveValue(page, el("Template input"), session.text("BDD Editor New {time}")));
      await session.step(92, "And \"Owner\" input should be visible", () => shouldBe(page, el("\"Owner\" input"), "visible"));
      await session.step(93, "And the \"Hit Design\" template \"BDD Editor New {time}\" on the server should have the key \"BDDEDN{time}\"", () => templateHas(page, "Hit Design", session.text("BDD Editor New {time}"), "key", session.text("BDDEDN{time}")));
      await session.step(94, "And the \"Hit Design\" template \"BDD Editor New {time}\" on the server should have the fields \"Owner (String, required), Due (Date)\"", () => templateHas(page, "Hit Design", session.text("BDD Editor New {time}"), "fields", "Owner (String, required), Due (Date)"));
      await session.step(95, "And the \"Hit Design\" template \"BDD Editor New {time}\" on the server should have the stages \"Design, Make\"", () => templateHas(page, "Hit Design", session.text("BDD Editor New {time}"), "stages", "Design, Make"));
      await session.step(96, "And the \"Hit Design\" template \"BDD Editor New {time}\" on the server should have the submit \"HitTriage:demoFileSubmit\"", () => templateHas(page, "Hit Design", session.text("BDD Editor New {time}"), "submit", "HitTriage:demoFileSubmit"));
      await session.step(97, "And the \"Hit Design\" template \"BDD Editor New {time}\" on the server should have the layout \"a layout\"", () => templateHas(page, "Hit Design", session.text("BDD Editor New {time}"), "layout", "a layout"));
      await session.step(98, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("CANCEL closes the editor and writes no template", async () => {
      await session.step(101, "When user closes all views", () => closeAllViews(page));
      await session.step(102, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(103, "And user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(104, "And user types \"BDD Editor Cancelled {time}\" into template name input", () => typeInto(page, session.text("BDD Editor Cancelled {time}"), el("template name input")));
      await session.step(105, "And user types \"BDDEDC{time}\" into template key input", () => typeInto(page, session.text("BDDEDC{time}"), el("template key input")));
      await session.step(106, "And user clicks on Cancel button", () => clickOn(page, el("Cancel button")));
      await session.step(107, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(108, "And Template input should offer the choice \"BDD Editor Existing {time}\"", () => shouldOfferChoice(page, el("Template input"), session.text("BDD Editor Existing {time}")));
      await session.step(109, "And Template input should not offer the choice \"BDD Editor Cancelled {time}\"", () => shouldNotOfferChoice(page, el("Template input"), session.text("BDD Editor Cancelled {time}")));
      await session.step(110, "And the \"Hit Design\" template \"BDD Editor Cancelled {time}\" should not be on the server", () => templateNotOnServer(page, "Hit Design", session.text("BDD Editor Cancelled {time}")));
    });
    await run.scenario("The crumb back to Hit Design closes the editor, and asks first once something was typed", async () => {
      await session.step(113, "When user closes all views", () => closeAllViews(page));
      await session.step(114, "And user opens the \"Hit Design\" app", () => openApp(page, "Hit Design"));
      await session.step(115, "And user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(116, "And user clicks on \"Hit Design\" link", () => clickOn(page, el("\"Hit Design\" link")));
      await session.step(117, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(118, "And the \"New Template\" view should not be open", () => viewNotOpen(page, "New Template"));
      await session.step(119, "When user clicks on \"New Template\" icon", () => clickOn(page, el("\"New Template\" icon")));
      await session.step(120, "And user types \"BDD Editor Cancelled {time}\" into template name input", () => typeInto(page, session.text("BDD Editor Cancelled {time}"), el("template name input")));
      await session.step(121, "And user clicks on \"Hit Design\" link", () => clickOn(page, el("\"Hit Design\" link")));
      await session.step(122, "Then \"Discard the new template?\" dialog should be visible", () => shouldBe(page, el("\"Discard the new template?\" dialog"), "visible"));
      await session.step(123, "When user clicks on Cancel button in \"Discard the new template?\" dialog", () => clickOn(page, el("Cancel button in \"Discard the new template?\" dialog")));
      await session.step(124, "Then the \"New Template\" view should be current", () => viewIsCurrent(page, "New Template"));
      await session.step(125, "And template name input should have the value \"BDD Editor Cancelled {time}\"", () => shouldHaveValue(page, el("template name input"), session.text("BDD Editor Cancelled {time}")));
      await session.step(126, "When user clicks on \"Hit Design\" link", () => clickOn(page, el("\"Hit Design\" link")));
      await session.step(127, "And user clicks on Discard button in \"Discard the new template?\" dialog", () => clickOn(page, el("Discard button in \"Discard the new template?\" dialog")));
      await session.step(128, "Then the \"Hit Design\" view should be current", () => viewIsCurrent(page, "Hit Design"));
      await session.step(129, "And the \"New Template\" view should not be open", () => viewNotOpen(page, "New Template"));
      await session.step(130, "And the \"Hit Design\" template \"BDD Editor Cancelled {time}\" should not be on the server", () => templateNotOnServer(page, "Hit Design", session.text("BDD Editor Cancelled {time}")));
    });
    run.finish();
  });
});
