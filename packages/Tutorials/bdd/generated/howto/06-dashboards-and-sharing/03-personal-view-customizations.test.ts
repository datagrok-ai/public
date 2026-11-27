/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/03-personal-view-customizations.feature
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
import {clickOn, close, enterInto, expand, hoverOver, selectIn, shouldBe, shouldContainText, shouldNotContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {balloonsClosed, browsePanelClosed, clickViewTab, contextPanelClosed, contextPanelShows, dialogCloses, doubleClickOpensProject, noProjectOnServer, openDataset, simpleModeOff, toolboxPaneHidden, userSettingsPutBack, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly, projectNotOpen} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {addViewerWith, noErrors, propertyShouldBe, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/03-personal-view-customizations.feature", import.meta.url);
  test("Change a dashboard just for yourself", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And simple mode is off", () => simpleModeOff(page));
    await session.step(24, "And the \"projectCustomViews\" user settings are put back at feature end", () => userSettingsPutBack(page, "projectCustomViews"));
    await session.step(25, "And the \"projectCustomLayout\" user settings are put back at feature end", () => userSettingsPutBack(page, "projectCustomLayout"));
    await session.step(26, "And no project named \"Demographics overview\" is on the server", () => noProjectOnServer(page, "Demographics overview"));
    await session.step(27, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(28, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(32, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["splitColumnName","RACE"]]), [["splitColumnName","RACE"]]);
    await session.step(34, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(35, "And user enters \"Demographics overview\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Demographics overview", el("Name text input in \"Save project\" dialog")));
    await session.step(36, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(37, "And the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(38, "And \"Share Demographics overview\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics overview\" dialog"), "visible"));
    await session.step(39, "And user clicks on CANCEL button in \"Share Demographics overview\" dialog", () => clickOn(page, el("CANCEL button in \"Share Demographics overview\" dialog")));
    await session.step(40, "And the \"Share Demographics overview\" dialog should close", () => dialogCloses(page, "Share Demographics overview"));
    await session.step(41, "And the notification balloons are closed", () => balloonsClosed(page));
    await session.step(42, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(43, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(44, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(45, "When the video says \"A team dashboard. You want it a bit different, just for you\"", () => videoSays(page, "A team dashboard. You want it a bit different, just for you"));
    await session.step(47, "And user picks \"RACE\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "RACE", "color", el("scatter plot viewer")), undefined, "Color the points by race instead of sex");
    await session.step(48, "Then \"colorColumnName\" property of scatter plot viewer should be \"RACE\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "RACE"));
    await session.step(49, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(50, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(52, "When user hovers over \"Save personal view customizations\" text in \"Save project\" dialog", () => hoverOver(page, el("\"Save personal view customizations\" text in \"Save project\" dialog")), undefined, "Save personal view customizations: for you only");
    await session.step(53, "And user selects \"Save personal view customizations\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save personal view customizations", el("radio input in \"Save project\" dialog")));
    await session.step(54, "Then Name text input in \"Save project\" dialog should be disabled", () => shouldBe(page, el("Name text input in \"Save project\" dialog"), "disabled"));
    await session.step(56, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")), undefined, "The dashboard itself stays as it was, for everyone else");
    await session.step(57, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(59, "When user closes \"demog\" view tab", () => close(page, el("\"demog\" view tab")), undefined, "Close it, and open it again later from Browse > Dashboards");
    await session.step(60, "Then the project \"Demographics overview\" should not be open", () => projectNotOpen(page, "Demographics overview"));
    await session.step(61, "And browse panel should be visible", () => shouldBe(page, el("browse panel"), "visible"));
    await session.step(62, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(63, "And user enters \"Demographics overview\" into gallery search", () => enterInto(page, "Demographics overview", el("gallery search")));
    await session.step(64, "Then \"Demographics overview\" gallery card should be visible", () => shouldBe(page, el("\"Demographics overview\" gallery card"), "visible"));
    await session.step(66, "When user double-clicks on \"Demographics overview\" gallery card and the project opens", () => doubleClickOpensProject(page, el("\"Demographics overview\" gallery card")), undefined, "It opens with your changes: colored by race");
    await session.step(67, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(68, "And a warning balloon containing \"Project has personal view customizations\" should have been shown", () => warningBalloonText(page, "Project has personal view customizations"));
    await session.step(69, "And the open table views should be exactly \"demog\"", () => openTableViewsExactly(page, "demog"));
    await session.step(70, "And \"colorColumnName\" property of scatter plot viewer should be \"RACE\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "RACE"));
    await session.step(72, "When user clicks on the tab of the \"Projects\" view", () => clickViewTab(page, "Projects"), undefined, "To go back to the shared version, select it in Dashboards");
    await session.step(73, "And user clicks on \"Demographics overview\" gallery card", () => clickOn(page, el("\"Demographics overview\" gallery card")));
    await session.step(75, "And user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "Its details are in the Context Panel");
    await session.step(76, "Then the context panel should show \"Demographics overview\"", () => contextPanelShows(page, "Demographics overview"));
    await session.step(78, "When user expands \"Custom views\" pane in context panel", () => expand(page, el("\"Custom views\" pane in context panel")), undefined, "Custom views: your version of this dashboard");
    await session.step(79, "Then \"Custom views\" pane in context panel should contain text \"Project has personal view customizations\"", () => shouldContainText(page, el("\"Custom views\" pane in context panel"), "Project has personal view customizations"));
    await session.step(81, "When user clicks on Reset button in \"Custom views\" pane in context panel", () => clickOn(page, el("Reset button in \"Custom views\" pane in context panel")), undefined, "Reset drops your version");
    await session.step(82, "Then \"Custom views\" pane in context panel should not contain text \"Project has personal view customizations\"", () => shouldNotContainText(page, el("\"Custom views\" pane in context panel"), "Project has personal view customizations"));
    await session.step(84, "When user closes \"demog\" view tab", () => close(page, el("\"demog\" view tab")), undefined, "Close it and open it again: the shared version, colored by sex");
    await session.step(85, "Then the project \"Demographics overview\" should not be open", () => projectNotOpen(page, "Demographics overview"));
    await session.step(86, "When user double-clicks on \"Demographics overview\" gallery card and the project opens", () => doubleClickOpensProject(page, el("\"Demographics overview\" gallery card")));
    await session.step(87, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(88, "And \"colorColumnName\" property of scatter plot viewer should be \"SEX\"", () => propertyShouldBe(page, "colorColumnName", el("scatter plot viewer"), "SEX"));
    await session.step(89, "And no errors should have been logged", () => noErrors(page));
  });
});
