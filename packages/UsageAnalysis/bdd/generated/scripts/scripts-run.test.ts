/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-run.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
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
import {clearField, clickOn, doubleClickOn, expand, selectIn, shouldBe, shouldHaveValue, typeInto, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {consoleCall, consoleShows, consoleShowsTimes, contextPanelOpen, contextPanelShows, dialogCloses, noteConsole, openDataset, scriptOnServer, scriptsView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Running a script with data from every source", () => {
  const session = feature(test, "features/scripts/scripts-run.feature", import.meta.url);
  test("Running a script with data from every source", {tag: ["@journey", "@full-stand", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And a script \"BddScriptRun{time}\" is on the server:", () => scriptOnServer(page, session.text("BddScriptRun{time}"), "#language: r\n#input: dataframe table\n#output: int count\n#output: string newParam\ncount <- nrow(table) * ncol(table)\nnewParam <- \"test\""));
    await session.step(32, "And the context panel is open", () => contextPanelOpen(page));
    await run.scenario("Run... with the open cars table counts its 510 cells", async () => {
      await session.step(35, "Given user opens cars dataset", () => openDataset(page, ds("cars")));
      await session.step(36, "And user opens the Scripts view", () => scriptsView(page));
      await session.step(37, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(38, "And user types \"BddScriptRun{time}\" into gallery search", () => typeInto(page, session.text("BddScriptRun{time}"), el("gallery search")));
      await session.step(39, "And user clicks on \"BddScriptRun{time}\" link in gallery", () => clickOn(page, el(session.text("\"BddScriptRun{time}\" link in gallery"))));
      await session.step(40, "Then the context panel should show \"BddScriptRun{time}\"", () => contextPanelShows(page, session.text("BddScriptRun{time}")));
      await session.step(41, "When user notes the console output", () => noteConsole(page));
      await session.step(42, "And user picks \"Run...\" from the context menu of \"BddScriptRun{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptRun{time}\" link in gallery"))));
      await session.step(43, "Then \"BddScriptRun{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BddScriptRun{time}\" dialog")), "visible"));
      await session.step(44, "When user selects \"cars\" in Table input in \"BddScriptRun{time}\" dialog", () => selectIn(page, "cars", el(session.text("Table input in \"BddScriptRun{time}\" dialog"))));
      await session.step(45, "And user clicks on OK button in \"BddScriptRun{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BddScriptRun{time}\" dialog"))));
      await session.step(46, "Then the \"BddScriptRun{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptRun{time}")));
      await session.step(47, "And the console should show \"count: 510\"", () => consoleShows(page, "count: 510"));
      await session.step(48, "And the console should show \"newParam: \\\"test\\\"\"", () => consoleShows(page, "newParam: \"test\""));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
      await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A local file is the table", async () => {
      await session.step(53, "When user notes the console output", () => noteConsole(page));
      await session.step(54, "And user picks \"Run...\" from the context menu of \"BddScriptRun{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptRun{time}\" link in gallery"))));
      await session.step(55, "Then \"BddScriptRun{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BddScriptRun{time}\" dialog")), "visible"));
      await session.step(56, "When user uploads \"fixtures/cars-small.csv\" through \"Open file\" icon in \"BddScriptRun{time}\" dialog", () => uploadThrough(page, "fixtures/cars-small.csv", el(session.text("\"Open file\" icon in \"BddScriptRun{time}\" dialog"))));
      await session.step(57, "Then Table input in \"BddScriptRun{time}\" dialog should have value \"cars-small\"", () => shouldHaveValue(page, el(session.text("Table input in \"BddScriptRun{time}\" dialog")), "cars-small"));
      await session.step(58, "When user clicks on OK button in \"BddScriptRun{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BddScriptRun{time}\" dialog"))));
      await session.step(59, "Then the \"BddScriptRun{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptRun{time}")));
      await session.step(60, "And the console should show \"count: 18\"", () => consoleShows(page, "count: 18"));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
      await session.step(62, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A file from Datagrok Files is the table", async () => {
      await session.step(65, "When user notes the console output", () => noteConsole(page));
      await session.step(66, "And user picks \"Run...\" from the context menu of \"BddScriptRun{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptRun{time}\" link in gallery"))));
      await session.step(67, "Then \"BddScriptRun{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BddScriptRun{time}\" dialog")), "visible"));
      await session.step(68, "When user clicks on \"Add file from Files\" icon in \"BddScriptRun{time}\" dialog", () => clickOn(page, el(session.text("\"Add file from Files\" icon in \"BddScriptRun{time}\" dialog"))));
      await session.step(69, "Then \"Select a file\" dialog should be visible", () => shouldBe(page, el("\"Select a file\" dialog"), "visible"));
      await session.step(70, "When user expands \"Files > Demo\" tree node inside \"Select a file\" dialog", () => expand(page, el("\"Files > Demo\" tree node inside \"Select a file\" dialog")));
      await session.step(71, "And user clicks on \"Files > Demo > cars.csv\" tree node inside \"Select a file\" dialog", () => clickOn(page, el("\"Files > Demo > cars.csv\" tree node inside \"Select a file\" dialog")));
      await session.step(72, "And user clicks on OK button in \"Select a file\" dialog", () => clickOn(page, el("OK button in \"Select a file\" dialog")));
      await session.step(73, "Then the \"Select a file\" dialog should close", () => dialogCloses(page, "Select a file"));
      await session.step(74, "And Table input in \"BddScriptRun{time}\" dialog should have value \"cars\"", () => shouldHaveValue(page, el(session.text("Table input in \"BddScriptRun{time}\" dialog")), "cars"));
      await session.step(75, "When user clicks on OK button in \"BddScriptRun{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BddScriptRun{time}\" dialog"))));
      await session.step(76, "Then the \"BddScriptRun{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptRun{time}")));
      await session.step(77, "And the console should show \"count: 510\"", () => consoleShows(page, "count: 510"));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
      await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A database query's result is the table", async () => {
      await session.step(82, "When user notes the console output", () => noteConsole(page));
      await session.step(83, "And user picks \"Run...\" from the context menu of \"BddScriptRun{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptRun{time}\" link in gallery"))));
      await session.step(84, "Then \"BddScriptRun{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BddScriptRun{time}\" dialog")), "visible"));
      await session.step(85, "When user clicks on \"Query database\" icon in \"BddScriptRun{time}\" dialog", () => clickOn(page, el(session.text("\"Query database\" icon in \"BddScriptRun{time}\" dialog"))));
      await session.step(86, "Then \"Select a database query\" dialog should be visible", () => shouldBe(page, el("\"Select a database query\" dialog"), "visible"));
      await session.step(87, "When user expands \"Postgres\" tree node inside \"Select a database query\" dialog", () => expand(page, el("\"Postgres\" tree node inside \"Select a database query\" dialog")));
      await session.step(88, "And user expands \"Postgres > NorthwindTest\" tree node inside \"Select a database query\" dialog", () => expand(page, el("\"Postgres > NorthwindTest\" tree node inside \"Select a database query\" dialog")));
      await session.step(89, "And user double-clicks on \"Postgres > NorthwindTest > PostgresAll\" tree node inside \"Select a database query\" dialog", () => doubleClickOn(page, el("\"Postgres > NorthwindTest > PostgresAll\" tree node inside \"Select a database query\" dialog")));
      await session.step(90, "Then the \"Select a database query\" dialog should close", () => dialogCloses(page, "Select a database query"));
      await session.step(91, "And Table input in \"BddScriptRun{time}\" dialog should have value \"PostgresAll\"", () => shouldHaveValue(page, el(session.text("Table input in \"BddScriptRun{time}\" dialog")), "PostgresAll"));
      await session.step(92, "When user clicks on OK button in \"BddScriptRun{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BddScriptRun{time}\" dialog"))));
      await session.step(93, "Then the \"BddScriptRun{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptRun{time}")));
      await session.step(94, "And the console should show \"count: 11620\"", () => consoleShows(page, "count: 11620"));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
      await session.step(96, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The console runs the script by its qualified name", async () => {
      await session.step(99, "When user notes the console output", () => noteConsole(page));
      await session.step(100, "And user calls the script \"BddScriptRun{time}\" from the console with '\"cars\"'", () => consoleCall(page, session.text("BddScriptRun{time}"), "\"cars\""));
      await session.step(101, "Then the console should show \"count: 510\"", () => consoleShows(page, "count: 510"));
      await session.step(102, "And the console should show \"newParam: \\\"test\\\"\"", () => consoleShows(page, "newParam: \"test\""));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("CANCEL runs nothing", async () => {
      await session.step(106, "Given user opens the Scripts view", () => scriptsView(page));
      await session.step(107, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(108, "And user types \"BddScriptRun{time}\" into gallery search", () => typeInto(page, session.text("BddScriptRun{time}"), el("gallery search")));
      await session.step(109, "And user notes the console output", () => noteConsole(page));
      await session.step(110, "And user picks \"Run...\" from the context menu of \"BddScriptRun{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptRun{time}\" link in gallery"))));
      await session.step(111, "Then \"BddScriptRun{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BddScriptRun{time}\" dialog")), "visible"));
      await session.step(112, "When user clicks on CANCEL button in \"BddScriptRun{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BddScriptRun{time}\" dialog"))));
      await session.step(113, "Then the \"BddScriptRun{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptRun{time}")));
      await session.step(115, "When user picks \"Run...\" from the context menu of \"BddScriptRun{time}\" link in gallery", () => pickFromContextMenu(page, "Run...", el(session.text("\"BddScriptRun{time}\" link in gallery"))));
      await session.step(116, "And user selects \"cars\" in Table input in \"BddScriptRun{time}\" dialog", () => selectIn(page, "cars", el(session.text("Table input in \"BddScriptRun{time}\" dialog"))));
      await session.step(117, "And user clicks on OK button in \"BddScriptRun{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BddScriptRun{time}\" dialog"))));
      await session.step(118, "Then the \"BddScriptRun{time}\" dialog should close", () => dialogCloses(page, session.text("BddScriptRun{time}")));
      await session.step(119, "And the console should show \"count: 510\" 1 time", () => consoleShowsTimes(page, "count: 510", 1));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
      await session.step(122, "When user clears gallery search", () => clearField(page, el("gallery search")));
    });
    run.finish();
  });
});
