/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/02-save-to-the-original-or-as-a-copy.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, close, enterInto, hoverOver, selectIn, shouldBe, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, closeAllViews, contextPanelClosed, dialogCloses, doubleClickOpensProject, noProjectOnServer, openDataset, openProjectWithTable, saveAsProject, simpleModeOff, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {projectNotOpen} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {addViewerWith, infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/02-save-to-the-original-or-as-a-copy.feature", import.meta.url);
  test("Save changes to the original or as a copy", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And no project named \"Demographics overview\" is on the server", () => noProjectOnServer(page, "Demographics overview"));
    await session.step(19, "And no project named \"Demographics by race\" is on the server", () => noProjectOnServer(page, "Demographics by race"));
    await session.step(20, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(21, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(25, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(27, "And user saves the current view as project \"Demographics overview\"", () => saveAsProject(page, "Demographics overview"));
    await session.step(28, "And user closes all views", () => closeAllViews(page));
    await session.step(29, "And user opens the \"Demographics overview\" project and waits for its table", () => openProjectWithTable(page, "Demographics overview"));
    await session.step(30, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(31, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(32, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(33, "When the video says \"This dashboard is already saved. Now change it\"", () => videoSays(page, "This dashboard is already saved. Now change it"));
    await session.step(35, "And user clicks on \"histogram\" icon in toolbox", () => clickOn(page, el("\"histogram\" icon in toolbox")), undefined, "Add a histogram");
    await session.step(36, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(38, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")), undefined, "SAVE now asks where the change should go");
    await session.step(39, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(41, "When user hovers over \"Save original project\" text in \"Save project\" dialog", () => hoverOver(page, el("\"Save original project\" text in \"Save project\" dialog")), undefined, "Save original project: everyone who opens it gets the change");
    await session.step(42, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(43, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(44, "And an info balloon containing 'Project \"Demographics overview\" uploaded' should have been shown", () => infoBalloonText(page, "Project \"Demographics overview\" uploaded"));
    await session.step(45, "When the video says \"For a variant, save a copy instead\"", () => videoSays(page, "For a variant, save a copy instead"));
    await session.step(47, "And user clicks on \"box plot\" icon in toolbox", () => clickOn(page, el("\"box plot\" icon in toolbox")), undefined, "Add a box plot");
    await session.step(48, "Then box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
    await session.step(49, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(50, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(52, "When user hovers over \"Save a copy\" text in \"Save project\" dialog", () => hoverOver(page, el("\"Save a copy\" text in \"Save project\" dialog")), undefined, "Save a copy: a new dashboard, the original stays as it was");
    await session.step(53, "And user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
    await session.step(54, "Then Name text input in \"Save project\" dialog should have value \"Copy of Demographics overview\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), "Copy of Demographics overview"));
    await session.step(56, "When user enters \"Demographics by race\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Demographics by race", el("Name text input in \"Save project\" dialog")), undefined, "Give the copy its own name");
    await session.step(57, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(58, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(59, "And an info balloon containing 'Project \"Demographics by race\" uploaded' should have been shown", () => infoBalloonText(page, "Project \"Demographics by race\" uploaded"));
    await session.step(61, "When user closes \"demog\" view tab", () => close(page, el("\"demog\" view tab")), undefined, "Close it. Browse > Dashboards now lists both");
    await session.step(62, "Then the project \"Demographics by race\" should not be open", () => projectNotOpen(page, "Demographics by race"));
    await session.step(63, "And browse panel should be visible", () => shouldBe(page, el("browse panel"), "visible"));
    await session.step(64, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(65, "And user enters \"Demographics\" into gallery search", () => enterInto(page, "Demographics", el("gallery search")));
    await session.step(66, "Then \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(67, "And \"Demographics by race\" gallery card should be visible", () => shouldBe(page, el("\"Demographics by race\" gallery card"), "visible"));
    await session.step(69, "When user hovers over \"Demographics by race\" gallery card", () => hoverOver(page, el("\"Demographics by race\" gallery card")), undefined, "The copy and the original, side by side");
    await session.step(71, "When user double-clicks on \"Demographics overview\" gallery card and the project opens", () => doubleClickOpensProject(page, el("\"Demographics overview\" gallery card")), undefined, "Open the original");
    await session.step(72, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(73, "And histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(74, "And box plot viewer should be absent", () => shouldBe(page, el("box plot viewer"), "absent"));
    await session.step(75, "When the video says \"It has the histogram, but not the box plot\"", () => videoSays(page, "It has the histogram, but not the box plot"));
    await session.step(76, "And user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
    await session.step(77, "Then no errors should have been logged", () => noErrors(page));
  });
});
