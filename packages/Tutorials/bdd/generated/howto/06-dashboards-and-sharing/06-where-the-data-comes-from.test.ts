/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/06-dashboards-and-sharing/06-where-the-data-comes-from.feature
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
import {clickOn, doubleClickOn, enterInto, hoverOver, isExpanded, shouldBe, shouldBeSwitchedOff, switchOff, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, browsePanelOpen, contextPanelClosed, dialogCloses, noProjectOnServer, savedAsSnapshot, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Dashboards", () => {
  const session = feature(test, "features/howto/06-dashboards-and-sharing/06-where-the-data-comes-from.feature", import.meta.url);
  test("Where a dashboard gets its data", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And no project named \"Demographics snapshot\" is on the server", () => noProjectOnServer(page, "Demographics snapshot"));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(20, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(21, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(22, "And user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(23, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(24, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["colorColumnName","SEX"]]);
    await session.step(28, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(29, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(30, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(31, "When the video says \"This table came from a file. A dashboard can read it again or keep a copy\"", () => videoSays(page, "This table came from a file. A dashboard can read it again or keep a copy"));
    await session.step(33, "And user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")), undefined, "SAVE lists the table with its Data sync switch");
    await session.step(34, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(35, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
    await session.step(37, "When user hovers over \"Data sync\" text in \"demog\" project table in \"Save project\" dialog", () => hoverOver(page, el("\"Data sync\" text in \"demog\" project table in \"Save project\" dialog")), undefined, "Data sync on: each time it opens, the dashboard reads the file again");
    await session.step(39, "And user clicks on \"Creation script\" button in \"demog\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog")), undefined, "Creation script: the steps that make the table again");
    await session.step(41, "And user switches off Data sync switch in \"demog\" project table in \"Save project\" dialog", () => switchOff(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog")), undefined, "Data sync off: the data is saved with the dashboard, as it is now");
    await session.step(42, "Then Data sync switch in \"demog\" project table in \"Save project\" dialog should be switched off", () => shouldBeSwitchedOff(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog")));
    await session.step(43, "When user enters \"Demographics snapshot\" into Name text input in \"Save project\" dialog", () => enterInto(page, "Demographics snapshot", el("Name text input in \"Save project\" dialog")));
    await session.step(44, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(45, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(46, "And the \"demog\" table of the \"Demographics snapshot\" project should be saved as a snapshot", () => savedAsSnapshot(page, "demog", "Demographics snapshot"));
    await session.step(47, "And \"Share Demographics snapshot\" dialog should be visible", () => shouldBe(page, el("\"Share Demographics snapshot\" dialog"), "visible"));
    await session.step(48, "When user clicks on CANCEL button in \"Share Demographics snapshot\" dialog", () => clickOn(page, el("CANCEL button in \"Share Demographics snapshot\" dialog")));
    await session.step(49, "Then the \"Share Demographics snapshot\" dialog should close", () => dialogCloses(page, "Share Demographics snapshot"));
    await session.step(50, "When the video says \"Keep it on for current data, off for a report that should not change\"", () => videoSays(page, "Keep it on for current data, off for a report that should not change"));
    await session.step(51, "And user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
    await session.step(52, "Then no errors should have been logged", () => noErrors(page));
  });
});
