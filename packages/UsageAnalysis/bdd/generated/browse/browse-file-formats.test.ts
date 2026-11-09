/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-file-formats.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
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
import {clickOn, doubleClickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewShows, packageInstalled, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The file formats of the demo files preview and open from Browse", () => {
  const session = feature(test, "features/browse/browse-file-formats.feature", import.meta.url);
  test("A CSV file previews as a table and opens into the workspace", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(29, "When user clicks on Files---Demo---cars.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
    await session.step(30, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(31, "And no table should be left in the workspace", () => noTableLeft(page));
    await session.step(32, "When user double-clicks on Files---Demo---cars.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
    await session.step(33, "Then table \"cars\" should have 30 rows", () => tableRows(page, "cars", 30));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
    await session.step(35, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A SDF file a package reads previews as a table and opens into the workspace [format=SDF, package=Chem, parent=Files---Demo---chem, folder=Files---Demo---chem---sdf, file=Files---Demo---chem---sdf---vitamin-d-2d.sdf, view=vitamin-d-2d, table=vitamin-d-2d, rows=1]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(38, "Given the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(39, "And Files---Demo---chem tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---chem tree node inside browse tree")));
    await session.step(40, "And Files---Demo---chem---sdf tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---chem---sdf tree node inside browse tree")));
    await session.step(41, "When user clicks on Files---Demo---chem---sdf---vitamin-d-2d.sdf tree node inside browse tree", () => clickOn(page, el("Files---Demo---chem---sdf---vitamin-d-2d.sdf tree node inside browse tree")));
    await session.step(42, "Then the \"vitamin-d-2d\" view should be current", () => viewIsCurrent(page, "vitamin-d-2d"));
    await session.step(43, "And no table should be left in the workspace", () => noTableLeft(page));
    await session.step(44, "When user double-clicks on Files---Demo---chem---sdf---vitamin-d-2d.sdf tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---chem---sdf---vitamin-d-2d.sdf tree node inside browse tree")));
    await session.step(45, "Then table \"vitamin-d-2d\" should have 1 rows", () => tableRows(page, "vitamin-d-2d", 1));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A Excel file a package reads previews as a table and opens into the workspace [format=Excel, package=PowerPack, parent=Files---Demo---northwind, folder=Files---Demo---northwind, file=Files---Demo---northwind---Northwind.xlsx, view=Suppliers, table=Orders, rows=830]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(38, "Given the \"PowerPack\" package is installed", () => packageInstalled(page, "PowerPack"));
    await session.step(39, "And Files---Demo---northwind tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---northwind tree node inside browse tree")));
    await session.step(40, "And Files---Demo---northwind tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---northwind tree node inside browse tree")));
    await session.step(41, "When user clicks on Files---Demo---northwind---Northwind.xlsx tree node inside browse tree", () => clickOn(page, el("Files---Demo---northwind---Northwind.xlsx tree node inside browse tree")));
    await session.step(42, "Then the \"Suppliers\" view should be current", () => viewIsCurrent(page, "Suppliers"));
    await session.step(43, "And no table should be left in the workspace", () => noTableLeft(page));
    await session.step(44, "When user double-clicks on Files---Demo---northwind---Northwind.xlsx tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---northwind---Northwind.xlsx tree node inside browse tree")));
    await session.step(45, "Then table \"Orders\" should have 830 rows", () => tableRows(page, "Orders", 830));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A ZIP file a package reads previews as a table and opens into the workspace [format=ZIP, package=PowerPack, parent=Files---Demo, folder=Files---Demo, file=Files---Demo---chika-mexicana.zip, view=Detailed Tickets from 2015-07-31 to 2015-08-07.xlsx - Worksheet, table=Detailed Tickets from 2015-05-08 to 2015-05-15.xlsx - Worksheet, rows=1859]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(38, "Given the \"PowerPack\" package is installed", () => packageInstalled(page, "PowerPack"));
    await session.step(39, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(40, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(41, "When user clicks on Files---Demo---chika-mexicana.zip tree node inside browse tree", () => clickOn(page, el("Files---Demo---chika-mexicana.zip tree node inside browse tree")));
    await session.step(42, "Then the \"Detailed Tickets from 2015-07-31 to 2015-08-07.xlsx - Worksheet\" view should be current", () => viewIsCurrent(page, "Detailed Tickets from 2015-07-31 to 2015-08-07.xlsx - Worksheet"));
    await session.step(43, "And no table should be left in the workspace", () => noTableLeft(page));
    await session.step(44, "When user double-clicks on Files---Demo---chika-mexicana.zip tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---chika-mexicana.zip tree node inside browse tree")));
    await session.step(45, "Then table \"Detailed Tickets from 2015-05-08 to 2015-05-15.xlsx - Worksheet\" should have 1859 rows", () => tableRows(page, "Detailed Tickets from 2015-05-08 to 2015-05-15.xlsx - Worksheet", 1859));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A JSON file previews as a document and stays open [format=JSON, parent=Files---Demo---kaggle, folder=Files---Demo---kaggle, file=Files---Demo---kaggle---insurance-metadata.json, view=insurance-metadata.json, text=@context]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(56, "Given Files---Demo---kaggle tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---kaggle tree node inside browse tree")));
    await session.step(57, "And Files---Demo---kaggle tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---kaggle tree node inside browse tree")));
    await session.step(58, "When user clicks on Files---Demo---kaggle---insurance-metadata.json tree node inside browse tree", () => clickOn(page, el("Files---Demo---kaggle---insurance-metadata.json tree node inside browse tree")));
    await session.step(59, "Then the \"insurance-metadata.json\" view should be current", () => viewIsCurrent(page, "insurance-metadata.json"));
    await session.step(60, "And the current view should show the text \"@context\"", () => currentViewShows(page, "@context"));
    await session.step(61, "When user double-clicks on Files---Demo---kaggle---insurance-metadata.json tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---kaggle---insurance-metadata.json tree node inside browse tree")));
    await session.step(62, "Then the \"insurance-metadata.json\" view should be current", () => viewIsCurrent(page, "insurance-metadata.json"));
    await session.step(63, "And the current view should show the text \"@context\"", () => currentViewShows(page, "@context"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
    await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A XML file previews as a document and stays open [format=XML, parent=Files---Demo---bio, folder=Files---Demo---bio---ngl-formats, file=Files---Demo---bio---ngl-formats---3pqr_validation.xml, view=3pqr_validation.xml, text=wwPDB-validation-information]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(56, "Given Files---Demo---bio tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---bio tree node inside browse tree")));
    await session.step(57, "And Files---Demo---bio---ngl-formats tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---bio---ngl-formats tree node inside browse tree")));
    await session.step(58, "When user clicks on Files---Demo---bio---ngl-formats---3pqr_validation.xml tree node inside browse tree", () => clickOn(page, el("Files---Demo---bio---ngl-formats---3pqr_validation.xml tree node inside browse tree")));
    await session.step(59, "Then the \"3pqr_validation.xml\" view should be current", () => viewIsCurrent(page, "3pqr_validation.xml"));
    await session.step(60, "And the current view should show the text \"wwPDB-validation-information\"", () => currentViewShows(page, "wwPDB-validation-information"));
    await session.step(61, "When user double-clicks on Files---Demo---bio---ngl-formats---3pqr_validation.xml tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---bio---ngl-formats---3pqr_validation.xml tree node inside browse tree")));
    await session.step(62, "Then the \"3pqr_validation.xml\" view should be current", () => viewIsCurrent(page, "3pqr_validation.xml"));
    await session.step(63, "And the current view should show the text \"wwPDB-validation-information\"", () => currentViewShows(page, "wwPDB-validation-information"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
    await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A HTML file previews as a document and stays open [format=HTML, parent=Files---Demo---forms, folder=Files---Demo---forms, file=Files---Demo---forms---tennis_players.html, view=tennis_players.html, text=display: flex]", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(56, "Given Files---Demo---forms tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---forms tree node inside browse tree")));
    await session.step(57, "And Files---Demo---forms tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---forms tree node inside browse tree")));
    await session.step(58, "When user clicks on Files---Demo---forms---tennis_players.html tree node inside browse tree", () => clickOn(page, el("Files---Demo---forms---tennis_players.html tree node inside browse tree")));
    await session.step(59, "Then the \"tennis_players.html\" view should be current", () => viewIsCurrent(page, "tennis_players.html"));
    await session.step(60, "And the current view should show the text \"display: flex\"", () => currentViewShows(page, "display: flex"));
    await session.step(61, "When user double-clicks on Files---Demo---forms---tennis_players.html tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---forms---tennis_players.html tree node inside browse tree")));
    await session.step(62, "Then the \"tennis_players.html\" view should be current", () => viewIsCurrent(page, "tennis_players.html"));
    await session.step(63, "And the current view should show the text \"display: flex\"", () => currentViewShows(page, "display: flex"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
    await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A Markdown file previews rendered and stays open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(74, "Given the \"PowerPack\" package is installed", () => packageInstalled(page, "PowerPack"));
    await session.step(75, "And Files---Demo---chem tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---chem tree node inside browse tree")));
    await session.step(76, "And Files---Demo---chem---protac tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---chem---protac tree node inside browse tree")));
    await session.step(77, "When user clicks on Files---Demo---chem---protac---protac_degraders.README.md tree node inside browse tree", () => clickOn(page, el("Files---Demo---chem---protac---protac_degraders.README.md tree node inside browse tree")));
    await session.step(78, "Then the \"protac_degraders\" view should be current", () => viewIsCurrent(page, "protac_degraders"));
    await session.step(79, "And the current view should show the text \"2,792 PROTAC degraders\"", () => currentViewShows(page, "2,792 PROTAC degraders"));
    await session.step(80, "When user double-clicks on Files---Demo---chem---protac---protac_degraders.README.md tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---chem---protac---protac_degraders.README.md tree node inside browse tree")));
    await session.step(81, "Then the \"protac_degraders\" view should be current", () => viewIsCurrent(page, "protac_degraders"));
    await session.step(82, "And the current view should show the text \"2,792 PROTAC degraders\"", () => currentViewShows(page, "2,792 PROTAC degraders"));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
    await session.step(84, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A PDF file previews its pages and stays open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(87, "Given the \"FileEditors\" package is installed", () => packageInstalled(page, "FileEditors"));
    await session.step(88, "And Files---Demo---texts tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---texts tree node inside browse tree")));
    await session.step(89, "And Files---Demo---texts---pdf tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---texts---pdf tree node inside browse tree")));
    await session.step(90, "When user clicks on Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree", () => clickOn(page, el("Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree")));
    await session.step(91, "Then the \"research-paper.pdf\" view should be current", () => viewIsCurrent(page, "research-paper.pdf"));
    await session.step(92, "And pdf preview page should be visible", () => shouldBe(page, el("pdf preview page"), "visible"));
    await session.step(93, "When user double-clicks on Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree")));
    await session.step(94, "Then the \"research-paper.pdf\" view should be current", () => viewIsCurrent(page, "research-paper.pdf"));
    await session.step(95, "And pdf preview page should be visible", () => shouldBe(page, el("pdf preview page"), "visible"));
    await session.step(96, "And no errors should have been logged", () => noErrors(page));
    await session.step(97, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A plate reader's text file previews as a plate with a tab per channel", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(25, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(26, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(100, "Given the \"Plates\" package is installed", () => packageInstalled(page, "Plates"));
    await session.step(101, "And Files---Demo---hts tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---hts tree node inside browse tree")));
    await session.step(102, "And Files---Demo---hts---formats tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo---hts---formats tree node inside browse tree")));
    await session.step(103, "When user clicks on Files---Demo---hts---formats---BMG_Pherastar_96.txt tree node inside browse tree", () => clickOn(page, el("Files---Demo---hts---formats---BMG_Pherastar_96.txt tree node inside browse tree")));
    await session.step(104, "Then the \"BMG_Pherastar_96.txt\" view should be current", () => viewIsCurrent(page, "BMG_Pherastar_96.txt"));
    await session.step(105, "And \"Channel A\" tab should be visible", () => shouldBe(page, el("\"Channel A\" tab"), "visible"));
    await session.step(106, "And no errors should have been logged", () => noErrors(page));
    await session.step(107, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
