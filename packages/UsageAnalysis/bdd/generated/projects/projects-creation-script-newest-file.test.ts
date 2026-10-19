/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-creation-script-newest-file.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
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
import {clickOn, enterInto, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewType, dialogCloses, doubleClickOpensProject, noProjectOnServer, reloadedByDataSync, scriptOnServer, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {copyInUsersFiles, noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A project whose table comes from the user's own script picks up the newest file on reopen", () => {
  const session = feature(test, "features/projects/projects-creation-script-newest-file.feature", import.meta.url);
  test("A project whose table comes from the user's own script picks up the newest file on reopen", {tag: ["@journey", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(19, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And no project named \"BDDCcsProject{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDCcsProject{time}")));
    await session.step(22, "And a copy of the \"System:DemoFiles/cars.csv\" file is in the user's files as \"BDDCcs{time}/data1.csv\"", () => copyInUsersFiles(page, "System:DemoFiles/cars.csv", session.text("BDDCcs{time}/data1.csv")));
    await session.step(23, "And a script \"BDDCcsScript{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDCcsScript{time}"), session.text("//language: javascript\n//output: dataframe result\nconst folder = grok.shell.user.project.name + ':Home/BDDCcs{time}/';\nconst csvFiles = await grok.dapi.files.list(folder, false, 'csv');\nif (csvFiles.length === 0)\n  throw new Error('No CSV files found in ' + folder);\nconst suffix = (name) => { const m = name.match(/(\\d+)(?=\\.csv$)/); return m ? parseInt(m[1], 10) : -1; };\ncsvFiles.sort((a, b) => suffix(a.fileName) - suffix(b.fileName));\nresult = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(csvFiles[csvFiles.length - 1].fullPath));")));
    await run.scenario("The script's result is saved with a bar chart, the script as its creation script", async () => {
      await session.step(37, "When user opens the address \"/func/Admin.BDDCcsScript{time}?run=true\"", () => openAddress(page, session.text("/func/Admin.BDDCcsScript{time}?run=true")));
      await session.step(38, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(39, "And the table should have 30 rows", () => rowCount(page, 30));
      await session.step(40, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(41, "When user clicks on \"bar chart\" icon in toolbox", () => clickOn(page, el("\"bar chart\" icon in toolbox")));
      await session.step(42, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
      await session.step(43, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(44, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(45, "And Data sync switch in \"BDDCcsScript{time}\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el(session.text("Data sync switch in \"BDDCcsScript{time}\" project table in \"Save project\" dialog")), "checked"));
      await session.step(46, "When user clicks on \"Creation script\" button in \"BDDCcsScript{time}\" project table in \"Save project\" dialog", () => clickOn(page, el(session.text("\"Creation script\" button in \"BDDCcsScript{time}\" project table in \"Save project\" dialog"))));
      await session.step(47, "Then \"BDDCcsScript{time}\" project table in \"Save project\" dialog should contain text \":BDDCcsScript{time}()\"", () => shouldContainText(page, el(session.text("\"BDDCcsScript{time}\" project table in \"Save project\" dialog")), session.text(":BDDCcsScript{time}()")));
      await session.step(48, "When user enters \"BDDCcsProject{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDCcsProject{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(49, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(50, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(51, "And an info balloon containing 'Project \"BDDCcsProject{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDCcsProject{time}\" uploaded")));
      await session.step(52, "And \"Share BDDCcsProject{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDCcsProject{time}\" dialog")), "visible"));
      await session.step(53, "When user clicks on CANCEL button in \"Share BDDCcsProject{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDCcsProject{time}\" dialog"))));
      await session.step(54, "Then the \"Share BDDCcsProject{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDCcsProject{time}")));
      await session.step(55, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(56, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(57, "And no table should be left in the workspace", () => noTableLeft(page));
    });
    await run.scenario("With a newer file in the folder, the reopened project shows its rows", async () => {
      await session.step(60, "Given a copy of the \"System:DemoFiles/iris.csv\" file is in the user's files as \"BDDCcs{time}/data2.csv\"", () => copyInUsersFiles(page, "System:DemoFiles/iris.csv", session.text("BDDCcs{time}/data2.csv")));
      await session.step(61, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(62, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(63, "And user enters \"BDDCcsProject{time}\" into gallery search", () => enterInto(page, session.text("BDDCcsProject{time}"), el("gallery search")));
      await session.step(64, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(65, "And user double-clicks on BDDCcsProject{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDCcsProject{time} gallery card"))));
      await session.step(66, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
      await session.step(67, "And the table should have 150 rows", () => rowCount(page, 150));
      await session.step(68, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(69, "And bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
      await session.step(70, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
      await session.step(71, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
