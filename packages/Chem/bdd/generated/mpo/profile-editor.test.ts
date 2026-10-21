/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/mpo/profile-editor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.mpo-profile-crud, chem.int.mpo-profile-sync]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {mpoProfileDescription, mpoProfileProperties, mpoProfilesOnServer, noMpoProfile, pmpoModelFile, shippedMpoProfile} from '../../bindings/mpo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, expand, isExpanded, selectIn, selectTextOf, shouldBe, shouldHaveValue, shouldOffer, typeInto, typeKeyByKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen, openDataset, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The MPO profile editor, from the Browse tree to a deleted profile", () => {
  const session = feature(test, "features/mpo/profile-editor.feature", import.meta.url);
  test("The MPO profile editor, from the Browse tree to a deleted profile", {tag: ["@serial", "@journey", "@realizes:chem.cp.mpo-profile-crud", "@realizes:chem.int.mpo-profile-sync"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(33, "And no MPO profile named \"BDD MPO {run}, BDD MPO {run} (Copy), BDD MPO renamed {run}, BDD MPO data-driven {run}\" is on the server", () => noMpoProfile(page, session.text("BDD MPO {run}, BDD MPO {run} (Copy), BDD MPO renamed {run}, BDD MPO data-driven {run}")));
    await session.step(34, "And the shipped MPO profile \"ADME Test\" is on the server", () => shippedMpoProfile(page, "ADME Test"));
    await session.step(35, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(37, "And Apps---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem tree node inside browse tree")));
    await run.scenario("The app opens from the Browse tree", async () => {
      await session.step(40, "When user clicks on Apps---Chem---MPO-profiles tree node inside browse tree", () => clickOn(page, el("Apps---Chem---MPO-profiles tree node inside browse tree")));
      await session.step(41, "Then the \"MPO Profiles\" view should be current", () => viewIsCurrent(page, "MPO Profiles"));
      await session.step(42, "And \"Create profile\" button should be visible", () => shouldBe(page, el("\"Create profile\" button"), "visible"));
      await session.step(43, "And 0 MPO profiles named \"BDD MPO {run}\" should be on the server", () => mpoProfilesOnServer(page, 0, session.text("BDD MPO {run}")));
      await session.step(44, "When user expands Apps---Chem---MPO-profiles tree node inside browse tree", () => expand(page, el("Apps---Chem---MPO-profiles tree node inside browse tree")));
      await session.step(45, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A property name typed key by key stays in its field", async () => {
      await session.step(48, "When user clicks on \"Create profile\" button", () => clickOn(page, el("\"Create profile\" button")));
      await session.step(49, "Then the \"Untitled Profile\" view should be current", () => viewIsCurrent(page, "Untitled Profile"));
      await session.step(50, "When user clicks on \"+ Add Property\" link", () => clickOn(page, el("\"+ Add Property\" link")));
      await session.step(51, "Then first MPO property should have value \"NewProperty1\"", () => shouldHaveValue(page, el("first MPO property"), "NewProperty1"));
      await session.step(52, "When user types \"HeavyAtomCount\" key by key into first MPO property", () => typeKeyByKey(page, "HeavyAtomCount", el("first MPO property")));
      await session.step(53, "Then first MPO property should have value \"HeavyAtomCount\"", () => shouldHaveValue(page, el("first MPO property"), "HeavyAtomCount"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Save stores the profile under the name typed on the tab", async () => {
      await session.step(57, "When user selects the text of MPO profile title", () => selectTextOf(page, el("MPO profile title")));
      await session.step(58, "And user types \"BDD MPO {run}\" into MPO profile title", () => typeInto(page, session.text("BDD MPO {run}"), el("MPO profile title")));
      await session.step(59, "And user clicks on \"Save\" button", () => clickOn(page, el("\"Save\" button")));
      await session.step(60, "Then an info balloon containing \"BDD MPO {run}\" should have been shown", () => infoBalloonText(page, session.text("BDD MPO {run}")));
      await session.step(61, "And 1 MPO profile named \"BDD MPO {run}\" should be on the server", () => mpoProfilesOnServer(page, 1, session.text("BDD MPO {run}")));
      await session.step(62, "And the MPO profile \"BDD MPO {run}\" should have the properties \"HeavyAtomCount\"", () => mpoProfileProperties(page, session.text("BDD MPO {run}"), "HeavyAtomCount"));
      await session.step(63, "And the \"BDD MPO {run}\" view should be current", () => viewIsCurrent(page, session.text("BDD MPO {run}")));
      await session.step(64, "And \"BDD MPO {run}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD MPO {run}\" tree node inside browse tree")), "visible"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The list holds the saved profile", async () => {
      await session.step(68, "When user switches to the \"MPO Profiles\" view", () => switchView(page, "MPO Profiles"));
      await session.step(69, "Then \"BDD MPO {run}\" MPO profile should be visible", () => shouldBe(page, el(session.text("\"BDD MPO {run}\" MPO profile")), "visible"));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A clone can be saved and then exists", async () => {
      await session.step(73, "When user clicks on actions of \"BDD MPO {run}\" MPO profile", () => clickOn(page, el(session.text("actions of \"BDD MPO {run}\" MPO profile"))));
      await session.step(74, "And user picks \"Clone\" from the open menu", () => pickFromOpenMenu(page, "Clone"));
      await session.step(75, "Then the \"BDD MPO {run} (Copy)\" view should be current", () => viewIsCurrent(page, session.text("BDD MPO {run} (Copy)")));
      await session.step(76, "And 0 MPO profiles named \"BDD MPO {run} (Copy)\" should be on the server", () => mpoProfilesOnServer(page, 0, session.text("BDD MPO {run} (Copy)")));
      await session.step(77, "And \"Save\" button should be enabled", () => shouldBe(page, el("\"Save\" button"), "enabled"));
      await session.step(78, "When user clicks on \"Save\" button", () => clickOn(page, el("\"Save\" button")));
      await session.step(79, "Then an info balloon containing \"BDD MPO {run} (Copy)\" should have been shown", () => infoBalloonText(page, session.text("BDD MPO {run} (Copy)")));
      await session.step(80, "And 1 MPO profile named \"BDD MPO {run} (Copy)\" should be on the server", () => mpoProfilesOnServer(page, 1, session.text("BDD MPO {run} (Copy)")));
      await session.step(81, "And the MPO profile \"BDD MPO {run} (Copy)\" should have the properties \"HeavyAtomCount\"", () => mpoProfileProperties(page, session.text("BDD MPO {run} (Copy)"), "HeavyAtomCount"));
      await session.step(82, "When user switches to the \"MPO Profiles\" view", () => switchView(page, "MPO Profiles"));
      await session.step(83, "Then \"BDD MPO {run} (Copy)\" MPO profile should be visible", () => shouldBe(page, el(session.text("\"BDD MPO {run} (Copy)\" MPO profile")), "visible"));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A profile renamed to a free name is saved without a replace prompt", async () => {
      await session.step(87, "When user clicks on actions of \"BDD MPO {run}\" MPO profile", () => clickOn(page, el(session.text("actions of \"BDD MPO {run}\" MPO profile"))));
      await session.step(88, "And user picks \"Edit\" from the open menu", () => pickFromOpenMenu(page, "Edit"));
      await session.step(89, "Then the \"BDD MPO {run}\" view should be current", () => viewIsCurrent(page, session.text("BDD MPO {run}")));
      await session.step(90, "When user selects the text of MPO profile title", () => selectTextOf(page, el("MPO profile title")));
      await session.step(91, "And user types \"BDD MPO renamed {run}\" into MPO profile title", () => typeInto(page, session.text("BDD MPO renamed {run}"), el("MPO profile title")));
      await session.step(92, "And user clicks on \"Save\" button", () => clickOn(page, el("\"Save\" button")));
      await session.step(93, "Then an info balloon containing \"BDD MPO renamed {run}\" should have been shown", () => infoBalloonText(page, session.text("BDD MPO renamed {run}")));
      await session.step(94, "And 1 MPO profile named \"BDD MPO renamed {run}\" should be on the server", () => mpoProfilesOnServer(page, 1, session.text("BDD MPO renamed {run}")));
      await session.step(95, "And 0 MPO profiles named \"BDD MPO {run}\" should be on the server", () => mpoProfilesOnServer(page, 0, session.text("BDD MPO {run}")));
      await session.step(96, "And the \"BDD MPO renamed {run}\" view should be current", () => viewIsCurrent(page, session.text("BDD MPO renamed {run}")));
      await session.step(97, "And \"BDD MPO renamed {run}\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD MPO renamed {run}\" tree node inside browse tree")), "visible"));
      await session.step(98, "And \"BDD MPO {run}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"BDD MPO {run}\" tree node inside browse tree")), "absent"));
      await session.step(99, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table opened after the profile tab is offered in its Dataset list", async () => {
      await session.step(102, "When user switches to the \"MPO Profiles\" view", () => switchView(page, "MPO Profiles"));
      await session.step(103, "And user clicks on \"Create profile\" button", () => clickOn(page, el("\"Create profile\" button")));
      await session.step(104, "Then Dataset input should have value \"\"", () => shouldHaveValue(page, el("Dataset input"), ""));
      await session.step(105, "Given user opens drugs-props-train dataset", () => openDataset(page, ds("drugs-props-train")));
      await session.step(106, "And user switches to the \"Untitled Profile\" view", () => switchView(page, "Untitled Profile"));
      await session.step(107, "Then Dataset input should offer \"drugs-props-train\"", () => shouldOffer(page, el("Dataset input"), "drugs-props-train"));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A data-driven profile is stored under the name and description typed on the tab", async () => {
      await session.step(111, "When user selects \"Data-driven\" in Method input", () => selectIn(page, "Data-driven", el("Method input")));
      await session.step(112, "And user selects \"drugs-props-train\" in Dataset input", () => selectIn(page, "drugs-props-train", el("Dataset input")));
      await session.step(113, "Then \"Training data\" heading should be visible", () => shouldBe(page, el("\"Training data\" heading"), "visible"));
      await session.step(114, "When user selects the text of MPO profile title", () => selectTextOf(page, el("MPO profile title")));
      await session.step(115, "And user types \"BDD MPO data-driven {run}\" into MPO profile title", () => typeInto(page, session.text("BDD MPO data-driven {run}"), el("MPO profile title")));
      await session.step(116, "And user types \"made from the CNS column\" into MPO profile description", () => typeInto(page, "made from the CNS column", el("MPO profile description")));
      await session.step(117, "And user clicks on \"Save\" button", () => clickOn(page, el("\"Save\" button")));
      await session.step(118, "Then an info balloon containing \"BDD MPO data-driven {run}\" should have been shown", () => infoBalloonText(page, session.text("BDD MPO data-driven {run}")));
      await session.step(119, "And 1 MPO profile named \"BDD MPO data-driven {run}\" should be on the server", () => mpoProfilesOnServer(page, 1, session.text("BDD MPO data-driven {run}")));
      await session.step(120, "And the MPO profile \"BDD MPO data-driven {run}\" should have the description \"made from the CNS column\"", () => mpoProfileDescription(page, session.text("BDD MPO data-driven {run}"), "made from the CNS column"));
      await session.step(121, "And the pMPO model file of \"BDD MPO data-driven {run}\" should hold its name and the description \"made from the CNS column\"", () => pmpoModelFile(page, session.text("BDD MPO data-driven {run}"), "made from the CNS column"));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Delete in the row menu removes the profile from the list and the server", async () => {
      await session.step(125, "When user switches to the \"MPO Profiles\" view", () => switchView(page, "MPO Profiles"));
      await session.step(126, "And user clicks on actions of \"BDD MPO renamed {run}\" MPO profile", () => clickOn(page, el(session.text("actions of \"BDD MPO renamed {run}\" MPO profile"))));
      await session.step(127, "And user picks \"Delete\" from the open menu", () => pickFromOpenMenu(page, "Delete"));
      await session.step(128, "Then \"Delete profile\" dialog should be visible", () => shouldBe(page, el("\"Delete profile\" dialog"), "visible"));
      await session.step(129, "When user clicks on OK button in \"Delete profile\" dialog", () => clickOn(page, el("OK button in \"Delete profile\" dialog")));
      await session.step(130, "Then 0 MPO profiles named \"BDD MPO renamed {run}\" should be on the server", () => mpoProfilesOnServer(page, 0, session.text("BDD MPO renamed {run}")));
      await session.step(131, "And \"BDD MPO renamed {run}\" MPO profile should be absent", () => shouldBe(page, el(session.text("\"BDD MPO renamed {run}\" MPO profile")), "absent"));
      await session.step(132, "And \"BDD MPO renamed {run}\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"BDD MPO renamed {run}\" tree node inside browse tree")), "absent"));
      await session.step(133, "And \"BDD MPO {run} (Copy)\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDD MPO {run} (Copy)\" tree node inside browse tree")), "visible"));
      await session.step(134, "When user clicks on actions of \"BDD MPO {run} (Copy)\" MPO profile", () => clickOn(page, el(session.text("actions of \"BDD MPO {run} (Copy)\" MPO profile"))));
      await session.step(135, "And user picks \"Delete\" from the open menu", () => pickFromOpenMenu(page, "Delete"));
      await session.step(136, "And user clicks on OK button in \"Delete profile\" dialog", () => clickOn(page, el("OK button in \"Delete profile\" dialog")));
      await session.step(137, "Then 0 MPO profiles named \"BDD MPO {run} (Copy)\" should be on the server", () => mpoProfilesOnServer(page, 0, session.text("BDD MPO {run} (Copy)")));
      await session.step(138, "And \"BDD MPO {run} (Copy)\" MPO profile should be absent", () => shouldBe(page, el(session.text("\"BDD MPO {run} (Copy)\" MPO profile")), "absent"));
      await session.step(139, "And \"BDD MPO {run} (Copy)\" tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("\"BDD MPO {run} (Copy)\" tree node inside browse tree")), "absent"));
      await session.step(140, "And \"ADME Test\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"ADME Test\" tree node inside browse tree"), "visible"));
      await session.step(141, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
