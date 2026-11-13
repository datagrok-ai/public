/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-dock-layout.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#239]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, dragTo, enterInto, isExpanded, pressKey, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, clickViewTab, closeAllViews, dialogCloses, dockViewBeside, doubleClickOpensProject, floatView, floatsAsRemembered, includeSaveDialogView, openDataset, openMarkupView, openProjectViews, openSaveDialog, ownProjectGone, projectsOnServer, rememberFloat, resizeViewShare, saveDialogViewExcluded, saveDialogViewIncluded, saveFromDialog, simpleModeOff, viewDockedBeside, viewFloating, viewTakesShare} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A dashboard saves the arrangement of its views", () => {
  const session = feature(test, "features/projects/projects-dock-layout.feature", import.meta.url);
  test("A dashboard saves the arrangement of its views", {tag: ["@journey", "@realizes:views.projects", "@realizes:core#239"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And simple mode is off", () => simpleModeOff(page));
    await session.step(21, "And the user's own project \"BDDLayout{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("BDDLayout{time}")));
    await run.scenario("Tables docked side by side and below, a markup view floating", async () => {
      await session.step(24, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(25, "And user opens iris dataset", () => openDataset(page, ds("iris")));
      await session.step(26, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(27, "When user drags the tab of the \"iris\" view to the right side of the \"cars\" view", () => dockViewBeside(page, "iris", "right", "cars"));
      await session.step(28, "Then the \"iris\" view should be docked right-of the \"cars\" view", () => viewDockedBeside(page, "iris", "right-of", "cars"));
      await session.step(29, "When user resizes the \"cars\" view to 30% of the width it shares with the \"iris\" view", () => resizeViewShare(page, "cars", 30, "width", "iris"));
      await session.step(30, "And user drags the tab of the \"demog-1000\" view to the bottom side of the \"cars\" view", () => dockViewBeside(page, "demog-1000", "bottom", "cars"));
      await session.step(31, "Then the \"demog-1000\" view should be docked below the \"cars\" view", () => viewDockedBeside(page, "demog-1000", "below", "cars"));
      await session.step(32, "And the \"iris\" view should be docked right-of the \"demog-1000\" view", () => viewDockedBeside(page, "iris", "right-of", "demog-1000"));
      await session.step(33, "And the \"cars\" view should take about 30% of the width it shares with the \"iris\" view", () => viewTakesShare(page, "cars", 30, "width", "iris"));
      await session.step(34, "Given user opens a markup view named \"BDD notes\"", () => openMarkupView(page, "BDD notes"));
      await session.step(35, "When user drags the tab of the \"BDD notes\" view out to float", () => floatView(page, "BDD notes"));
      await session.step(36, "Then the \"BDD notes\" view should be floating", () => viewFloating(page, "BDD notes"));
      await session.step(37, "When user remembers where the \"BDD notes\" view floats", () => rememberFloat(page, "BDD notes"));
      await session.step(38, "Given user opens a markup view named \"BDD scratch\"", () => openMarkupView(page, "BDD scratch"));
      await session.step(39, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Save dialog lists the markup views left out until one is included", async () => {
      await session.step(42, "When user clicks on the tab of the \"cars\" view", () => clickViewTab(page, "cars"));
      await session.step(43, "And user opens the Save project dialog from the ribbon", () => openSaveDialog(page));
      await session.step(44, "Then the \"BDD notes\" view should be excluded in the Save project dialog", () => saveDialogViewExcluded(page, "BDD notes"));
      await session.step(45, "And the \"BDD scratch\" view should be excluded in the Save project dialog", () => saveDialogViewExcluded(page, "BDD scratch"));
      await session.step(46, "When user includes the \"BDD notes\" view in the Save project dialog", () => includeSaveDialogView(page, "BDD notes"));
      await session.step(47, "Then the \"BDD scratch\" view should be excluded in the Save project dialog", () => saveDialogViewExcluded(page, "BDD scratch"));
      await session.step(48, "When user types \"BDDLayout{time}\" into text input in \"Save project\" dialog", () => typeInto(page, session.text("BDDLayout{time}"), el("text input in \"Save project\" dialog")));
      await session.step(49, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(50, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(51, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(52, "Then 1 project named \"BDDLayout{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDLayout{time}")));
    });
    await run.scenario("The reopened dashboard brings the arrangement back", async () => {
      await session.step(55, "When user closes all views", () => closeAllViews(page));
      await session.step(56, "Then no table should be left in the workspace", () => noTableLeft(page));
      await session.step(57, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(58, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(59, "And user enters \"BDDLayout{time}\" into gallery search", () => enterInto(page, session.text("BDDLayout{time}"), el("gallery search")));
      await session.step(60, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(61, "And user double-clicks on BDDLayout{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDLayout{time} gallery card"))));
      await session.step(62, "Then the \"iris\" view should be docked right-of the \"cars\" view", () => viewDockedBeside(page, "iris", "right-of", "cars"));
      await session.step(63, "And the \"cars\" view should take about 30% of the width it shares with the \"iris\" view", () => viewTakesShare(page, "cars", 30, "width", "iris"));
      await session.step(64, "And the \"demog-1000\" view should be docked below the \"cars\" view", () => viewDockedBeside(page, "demog-1000", "below", "cars"));
      await session.step(65, "And the \"iris\" view should be docked right-of the \"demog-1000\" view", () => viewDockedBeside(page, "iris", "right-of", "demog-1000"));
      await session.step(66, "And the \"BDD notes\" view should float where it was remembered", () => floatsAsRemembered(page, "BDD notes"));
      await session.step(67, "And \"BDD scratch\" view should be absent", () => shouldBe(page, el("\"BDD scratch\" view"), "absent"));
      await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A view dragged onto the open dashboard in the Dashboards panel is saved with it", async () => {
      await session.step(71, "Given user opens a markup view named \"BDD moved\"", () => openMarkupView(page, "BDD moved"));
      await session.step(72, "When user clicks on Dashboards tab", () => clickOn(page, el("Dashboards tab")));
      await session.step(73, "Then \"New Dashboard > BDD moved\" tree node inside browse tree should be visible", () => shouldBe(page, el("\"New Dashboard > BDD moved\" tree node inside browse tree"), "visible"));
      await session.step(74, "When user drags \"New Dashboard > BDD moved\" tree node inside browse tree to BDDLayout{time} tree node inside browse tree", () => dragTo(page, el("\"New Dashboard > BDD moved\" tree node inside browse tree"), el(session.text("BDDLayout{time} tree node inside browse tree"))));
      await session.step(75, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(76, "And Move entity dialog should contain text \"BDD moved\"", () => shouldContainText(page, el("Move entity dialog"), "BDD moved"));
      await session.step(77, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(78, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(79, "Given BDDLayout{time} tree node inside browse tree is expanded", () => isExpanded(page, el(session.text("BDDLayout{time} tree node inside browse tree"))));
      await session.step(80, "Then \"BDDLayout{time} > BDD moved\" tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("\"BDDLayout{time} > BDD moved\" tree node inside browse tree")), "visible"));
      await session.step(81, "When user clicks on Save button in BDDLayout{time} tree node inside browse tree", () => clickOn(page, el(session.text("Save button in BDDLayout{time} tree node inside browse tree"))));
      await session.step(82, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(83, "And the \"BDD moved\" view should be included in the Save project dialog", () => saveDialogViewIncluded(page, "BDD moved"));
      await session.step(84, "When user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
      await session.step(85, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(86, "When user closes all views", () => closeAllViews(page));
      await session.step(87, "Then no table should be left in the workspace", () => noTableLeft(page));
      await session.step(88, "When user opens the \"BDDLayout{time}\" project and waits until it has opened", () => openProjectViews(page, session.text("BDDLayout{time}")));
      await session.step(89, "Then \"BDD moved\" view should be present", () => shouldBe(page, el("\"BDD moved\" view"), "present"));
      await session.step(90, "And the \"iris\" view should be docked right-of the \"cars\" view", () => viewDockedBeside(page, "iris", "right-of", "cars"));
      await session.step(91, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
