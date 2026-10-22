/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/moltrack-ems-registration.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {compoundCount, hasRememberedText, namesRemembered, noCompound, openEms, rememberCompound, rememberText} from '../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldNotContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {shouldHaveItemsInTotal} from '@datagrok-libraries/bdd/bindings/common/u2-steps';
import {packageInstalled} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, enter, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("MolTrack EMS: registering a compound", () => {
  const session = feature(test, "features/moltrack-ems-registration.feature", import.meta.url);
  test("MolTrack EMS: registering a compound", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(12, "And no compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" is registered, now and at feature end", () => noCompound(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"));
    await session.step(13, "And user opens the MolTrack EMS app", () => openEms(page));
    enter(page, "MolTrack EMS");
    await run.scenario("A new structure is registered through the form", async () => {
      await session.step(16, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")));
      await session.step(17, "Then breadcrumbs should contain text \"New compound\"", () => shouldContainText(page, el("breadcrumbs"), "New compound"));
      await session.step(18, "When user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")));
      await session.step(19, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O", el("molecule input of sketcher dialog")));
      await session.step(20, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(21, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(22, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
      await session.step(23, "And \"Formula\" read-only field should contain text \"C28H56O2\"", () => shouldContainText(page, el("\"Formula\" read-only field"), "C28H56O2"));
      await session.step(24, "When user types \"Octacosanoic acid\" into \"Common Name\" text input", () => typeInto(page, "Octacosanoic acid", el("\"Common Name\" text input")));
      await session.step(25, "And user clicks on \"Save\" button in app ribbon", () => clickOn(page, el("\"Save\" button in app ribbon")));
      await session.step(26, "Then notification should contain text \"Compound saved\"", () => shouldContainText(page, el("notification"), "Compound saved"));
      await session.step(27, "And breadcrumbs should not contain text \"New compound\"", () => shouldNotContainText(page, el("breadcrumbs"), "New compound"));
      await session.step(28, "And 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"));
      await session.step(29, "When user remembers the compound the breadcrumbs name", () => rememberCompound(page));
    });
    await run.scenario("Its sodium salt is refused as the compound already registered", async () => {
      await session.step(32, "When user clicks on \"New\" button in app ribbon", () => clickOn(page, el("\"New\" button in app ribbon")));
      await session.step(33, "And user clicks on molecule canvas in \"Structure\" dart input", () => clickOn(page, el("molecule canvas in \"Structure\" dart input")));
      await session.step(34, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)[O-].[Na+]", el("molecule input of sketcher dialog")));
      await session.step(35, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(36, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(37, "Then \"Structure\" dart input should contain text \"This structure is already registered as DG-\"", () => shouldContainText(page, el("\"Structure\" dart input"), "This structure is already registered as DG-"));
      await session.step(38, "And \"Structure\" dart input should name the remembered compound", () => namesRemembered(page, el("\"Structure\" dart input")));
      await session.step(39, "And \"Save\" button in app ribbon should be disabled", () => shouldBe(page, el("\"Save\" button in app ribbon"), "disabled"));
      await session.step(40, "When user clicks on \"Discard\" button in app ribbon", () => clickOn(page, el("\"Discard\" button in app ribbon")));
      await session.step(41, "Then domain list should be visible", () => shouldBe(page, el("domain list"), "visible"));
      await session.step(42, "And 1 compound with the structure \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O\" should be registered on the server", () => compoundCount(page, 1, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)O"));
    });
    await run.scenario("A substructure condition narrows the list to the compound", async () => {
      await session.step(45, "When user remembers the text of compound count", () => rememberText(page, el("compound count")));
      await session.step(46, "And user clicks on \"Show the filters\" button", () => clickOn(page, el("\"Show the filters\" button")));
      await session.step(47, "And user clicks on add condition in filters pane", () => clickOn(page, el("add condition in filters pane")));
      await session.step(48, "And user selects \"Parent Structure\" in condition property", () => selectIn(page, "Parent Structure", el("condition property")));
      await session.step(49, "And user selects \"has substructure\" in condition operator", () => selectIn(page, "has substructure", el("condition operator")));
      await session.step(50, "And user clicks on molecule canvas in condition value", () => clickOn(page, el("molecule canvas in condition value")));
      await session.step(51, "And user types \"CCCCCCCCCCCCCCCCCCCCCCCCCCCC\" into molecule input of sketcher dialog", () => typeInto(page, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", el("molecule input of sketcher dialog")));
      await session.step(52, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(53, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(54, "Then compound list should have 1 item in total", () => shouldHaveItemsInTotal(page, el("compound list"), 1));
      await session.step(55, "And compound count should have text \"1 compound\"", () => shouldHaveText(page, el("compound count"), "1 compound"));
      await session.step(56, "And compound list should contain text \"Octacosanoic acid\"", () => shouldContainText(page, el("compound list"), "Octacosanoic acid"));
      await session.step(57, "And compound list should name the remembered compound", () => namesRemembered(page, el("compound list")));
      await session.step(58, "When user clicks on condition remove", () => clickOn(page, el("condition remove")));
      await session.step(59, "Then compound count should have the remembered text", () => hasRememberedText(page, el("compound count")));
    });
    run.finish();
  });
});
