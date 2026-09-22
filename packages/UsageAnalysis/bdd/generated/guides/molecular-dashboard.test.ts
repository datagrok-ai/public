/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/molecular-dashboard.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noProjectOnServer, openDataset, projectsOnServer, simpleModeOff, viewHoldsViewers} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, openToolbox, paintedInColors, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Build a dashboard on molecular data", () => {
  const session = feature(test, "features/guides/molecular-dashboard.feature", import.meta.url);
  test("Add viewers from the toolbox and set each one up", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And no project named \"Molecular dashboard\" is on the server", () => noProjectOnServer(page, "Molecular dashboard"));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "When user opens toolbox", () => openToolbox(page));
    await session.step(17, "And user clicks on scatter plot icon on toolbox", () => clickOn(page, el("scatter plot icon on toolbox")));
    await session.step(18, "Then scatter plot viewer should be added to the open tableview", () => viewerAdded(page, "scatter plot"));
    await session.step(19, "When user picks \"Chemical Space X\" in the \"x\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "Chemical Space X", "x", el("scatter plot viewer")));
    await session.step(20, "And user picks \"Chemical Space Y\" in the \"y\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "Chemical Space Y", "y", el("scatter plot viewer")));
    await session.step(21, "And user picks \"Series\" in the \"color\" column selector of scatter plot viewer", () => pickInColumnSelector(page, "Series", "color", el("scatter plot viewer")));
    await session.step(22, "Then scatter plot viewer should be painted in at least 3 colors", () => paintedInColors(page, el("scatter plot viewer"), 3));
    await session.step(23, "When user clicks on histogram icon on toolbox", () => clickOn(page, el("histogram icon on toolbox")));
    await session.step(24, "And user selects \"Average Mass\" in Value column input in histogram viewer", () => selectIn(page, "Average Mass", el("Value column input in histogram viewer")));
    await session.step(25, "And user clicks on bar chart icon on toolbox", () => clickOn(page, el("bar chart icon on toolbox")));
    await session.step(26, "And user picks \"Series\" in the \"split\" column selector of bar chart viewer", () => pickInColumnSelector(page, "Series", "split", el("bar chart viewer")));
    await session.step(27, "And user clicks on pie chart icon on toolbox", () => clickOn(page, el("pie chart icon on toolbox")));
    await session.step(28, "And user picks \"Stereo Category\" in the \"category\" column selector of pie chart viewer", () => pickInColumnSelector(page, "Stereo Category", "category", el("pie chart viewer")));
    await session.step(29, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(30, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
    await session.step(31, "When user clicks on settings icon of bar chart viewer", () => clickOn(page, el("settings icon of bar chart viewer")));
    await session.step(32, "Given \"Y Axis\" category in context panel is expanded", () => isExpanded(page, el("\"Y Axis\" category in context panel")));
    await session.step(33, "When user selects \"avg\" in \"Value Aggr Type\" property in context panel", () => selectIn(page, "avg", el("\"Value Aggr Type\" property in context panel")));
    await session.step(34, "And user selects \"Average Mass\" in \"Value\" property in context panel", () => selectIn(page, "Average Mass", el("\"Value\" property in context panel")));
    await session.step(35, "Then the current view should hold at least 5 viewers", () => viewHoldsViewers(page, 5));
    await session.step(36, "When user clicks on Save button in toolbar", () => clickOn(page, el("Save button in toolbar")));
    await session.step(37, "And user enters \"Molecular dashboard\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Molecular dashboard", el("Name text input in \"Save project\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(39, "Then 1 project named \"Molecular dashboard\" should be on the server", () => projectsOnServer(page, 1, "Molecular dashboard"));
    await session.step(40, "When user clicks on CANCEL button in \"Share Molecular dashboard\" dialog", () => clickOn(page, el("CANCEL button in \"Share Molecular dashboard\" dialog")));
    await session.step(41, "Then no errors should have been logged", () => noErrors(page));
  });
});
