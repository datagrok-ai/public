/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-creation-script.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noErrorsButPreviewNoise} from '../../bindings/projects-regressions.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe, shouldBeSwitchedOn, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, dialogCloses, noProjectOnServer, projectsOnServer, savedWithDataSync, scriptOnServer, scriptsView, switchTableView, tableReloaded, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {homeFileWritten, noHomeFolder, noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project whose table comes from a script picks up the newest file on reopen", () => {
  const session = feature(test, "features/projects/projects-creation-script.feature", import.meta.url);
  test("A project whose table comes from a script picks up the newest file on reopen", {tag: ["@journey", "@serial", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And no project named \"BDDCreate{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDCreate{time}")));
    await session.step(27, "And no folder \"BDDCreate{time}\" is in the user's files", () => noHomeFolder(page, session.text("BDDCreate{time}")));
    await session.step(28, "And the file \"data_1.csv\" of the folder \"BDDCreate{time}\" in the user's files holds:", () => homeFileWritten(page, "data_1.csv", session.text("BDDCreate{time}"), "source,value\ndata_1,1\ndata_1,2\ndata_1,3"));
    await session.step(35, "And a script \"BDDCreate{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDCreate{time}"), session.text("//language: javascript\n//output: dataframe result\nconst folder = grok.shell.user.project.name + ':Home/BDDCreate{time}/';\nconst csvFiles = await grok.dapi.files.list(folder, false, 'csv');\nif (csvFiles.length === 0)\n  throw new Error('No CSV files found in ' + folder);\nconst suffix = (name) => { const m = name.match(/(\\d+)(?=\\.csv$)/); return m ? parseInt(m[1], 10) : -1; };\ncsvFiles.sort((a, b) => suffix(a.fileName) - suffix(b.fileName));\nresult = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(csvFiles[csvFiles.length - 1].fullPath));")));
    await run.scenario("The script run from Browse returns the highest-suffixed file", async () => {
      await session.step(49, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
      await session.step(50, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
      await session.step(51, "When user clicks on Platform---Functions---Scripts tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Scripts tree node inside browse tree")));
      await session.step(52, "Then the \"Scripts\" view should be current", () => viewIsCurrent(page, "Scripts"));
      await session.step(53, "When user enters \"BDDCreate{time}\" into gallery search", () => enterInto(page, session.text("BDDCreate{time}"), el("gallery search")));
      await session.step(54, "And user picks \"Run...\" from the context menu of BDDCreate{time} gallery card", () => pickFromContextMenu(page, "Run...", el(session.text("BDDCreate{time} gallery card"))));
      await session.step(55, "Then the \"result\" view should be current", () => viewIsCurrent(page, "result"));
      await session.step(56, "And table \"result\" should have 3 rows", () => tableRows(page, "result", 3));
      await session.step(57, "And every value of \"source\" column should match \"^data_1$\"", () => everyValueMatches(page, "source", "^data_1$"));
      await session.step(59, "Given user opens the Scripts view", () => scriptsView(page));
      await session.step(60, "And user switches to the \"result\" table view", () => switchTableView(page, "result"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
      await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The result is saved with Data sync, the script as its creation script", async () => {
      await session.step(65, "When user clicks on \"bar chart\" icon in toolbox", () => clickOn(page, el("\"bar chart\" icon in toolbox")));
      await session.step(66, "Then the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
      await session.step(67, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(68, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(69, "And Data sync switch in \"result\" project table in \"Save project\" dialog should be switched on", () => shouldBeSwitchedOn(page, el("Data sync switch in \"result\" project table in \"Save project\" dialog")));
      await session.step(70, "When user clicks on \"Creation script\" button in \"result\" project table in \"Save project\" dialog", () => clickOn(page, el("\"Creation script\" button in \"result\" project table in \"Save project\" dialog")));
      await session.step(71, "Then \"result\" project table in \"Save project\" dialog should contain text \"BDDCreate{time}\"", () => shouldContainText(page, el("\"result\" project table in \"Save project\" dialog"), session.text("BDDCreate{time}")));
      await session.step(72, "When user enters \"BDDCreate{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDCreate{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(73, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(74, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(75, "And an info balloon containing 'Project \"BDDCreate{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDCreate{time}\" uploaded")));
      await session.step(76, "And 1 project named \"BDDCreate{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("BDDCreate{time}")));
      await session.step(77, "And the \"result\" table of the \"BDDCreate{time}\" project should be saved with data sync", () => savedWithDataSync(page, "result", session.text("BDDCreate{time}")));
      await session.step(78, "And \"Share BDDCreate{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDCreate{time}\" dialog")), "visible"));
      await session.step(79, "When user clicks on CANCEL button in \"Share BDDCreate{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDCreate{time}\" dialog"))));
      await session.step(80, "Then the \"Share BDDCreate{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDCreate{time}")));
      await session.step(81, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(82, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(83, "And no table should be left in the workspace", () => noTableLeft(page));
      await session.step(84, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
    });
    await run.scenario("A newer file in the folder is what the reopened project shows", async () => {
      await session.step(87, "Given the file \"data_2.csv\" of the folder \"BDDCreate{time}\" in the user's files holds:", () => homeFileWritten(page, "data_2.csv", session.text("BDDCreate{time}"), "source,value\ndata_2,10\ndata_2,20\ndata_2,30\ndata_2,40\ndata_2,50"));
      await session.step(96, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(97, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(98, "And user enters \"BDDCreate{time}\" into gallery search", () => enterInto(page, session.text("BDDCreate{time}"), el("gallery search")));
      await session.step(99, "And user double-clicks on BDDCreate{time} gallery card", () => doubleClickOn(page, el(session.text("BDDCreate{time} gallery card"))));
      await session.step(100, "Then table \"result\" should have been reloaded by data sync with 5 rows", () => tableReloaded(page, "result", 5));
      await session.step(101, "And every value of \"source\" column should match \"^data_2$\"", () => everyValueMatches(page, "source", "^data_2$"));
      await session.step(102, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
      await session.step(104, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
