/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/connections/connections-import-swagger.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, enterSecret, followingShouldBe, isExpanded, openLocalFile, rightClickOn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, closeAllViews, connectionsOnServer, dialogCloses, noConnectionOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noBalloons, noErrors, openContextMenu, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Importing an OpenAPI (Swagger) file as a connection", () => {
  const session = feature(test, "features/connections/connections-import-swagger.feature", import.meta.url);
  test("The imported file is a connection with a query per operation", {tag: ["@connections", "@full-stand", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And no connection named \"BDD-Conn-Swagger\" is on the server", () => noConnectionOnServer(page, "BDD-Conn-Swagger"));
    await session.step(27, "When user opens the local file \"fixtures/bdd-swagger.yaml\"", () => openLocalFile(page, "fixtures/bdd-swagger.yaml"));
    await session.step(28, "Then 1 connection named \"BDD-Conn-Swagger\" should be on the server", () => connectionsOnServer(page, 1, "BDD-Conn-Swagger"));
    await session.step(29, "And the \"bdd-swagger.yaml\" view should be current", () => viewIsCurrent(page, "bdd-swagger.yaml"));
    await session.step(30, "When user closes all views", () => closeAllViews(page));
    await session.step(31, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(34, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(35, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(36, "And Platform---Functions---OpenAPI tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions---OpenAPI tree node inside browse tree")));
    await session.step(39, "Then Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree should be visible", () => shouldBe(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree"), "visible"));
    await session.step(40, "Given Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree")));
    await session.step(41, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Current-Weather-Data-By-City-Name tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Cities-In-Cycle tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---5-day/3-hour-Forecast-By-City-Name tree node inside browse tree"]]), [["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Current-Weather-Data-By-City-Name tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Cities-In-Cycle tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---5-day/3-hour-Forecast-By-City-Name tree node inside browse tree"]]);
    await session.step(45, "When user opens the context menu of Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree", () => openContextMenu(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree")));
    await session.step(46, "Then the open menu should list \"Edit...\"", () => menuLists(page, "Edit..."));
    await session.step(47, "And the open menu should list \"Test connection\"", () => menuLists(page, "Test connection"));
    await session.step(48, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(49, "Then no errors should have been logged", () => noErrors(page));
  });
  test("With its API key the connection's query returns the weather", {tag: ["@connections", "@full-stand", "@serial", "@needs-credentials"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And no connection named \"BDD-Conn-Swagger\" is on the server", () => noConnectionOnServer(page, "BDD-Conn-Swagger"));
    await session.step(27, "When user opens the local file \"fixtures/bdd-swagger.yaml\"", () => openLocalFile(page, "fixtures/bdd-swagger.yaml"));
    await session.step(28, "Then 1 connection named \"BDD-Conn-Swagger\" should be on the server", () => connectionsOnServer(page, 1, "BDD-Conn-Swagger"));
    await session.step(29, "And the \"bdd-swagger.yaml\" view should be current", () => viewIsCurrent(page, "bdd-swagger.yaml"));
    await session.step(30, "When user closes all views", () => closeAllViews(page));
    await session.step(31, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(34, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(35, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(36, "And Platform---Functions---OpenAPI tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions---OpenAPI tree node inside browse tree")));
    await session.step(53, "When user right-clicks on Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree", () => rightClickOn(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree")));
    await session.step(54, "And user picks \"Edit...\" from the open menu", () => pickFromOpenMenu(page, "Edit..."));
    await session.step(55, "Then \"Edit Connection\" dialog should be visible", () => shouldBe(page, el("\"Edit Connection\" dialog"), "visible"));
    await session.step(56, "And Url input in \"Edit Connection\" dialog should have value \"https://api.openweathermap.org/data/2.5\"", () => shouldHaveValue(page, el("Url input in \"Edit Connection\" dialog"), "https://api.openweathermap.org/data/2.5"));
    await session.step(57, "When user enters the DG_OPENWEATHERMAP_API_KEY secret into ApiKey input in \"Edit Connection\" dialog", () => enterSecret(page, "DG_OPENWEATHERMAP_API_KEY", el("ApiKey input in \"Edit Connection\" dialog")));
    await session.step(58, "And user clicks on OK button in \"Edit Connection\" dialog", () => clickOn(page, el("OK button in \"Edit Connection\" dialog")));
    await session.step(59, "Then the \"Edit Connection\" dialog should close", () => dialogCloses(page, "Edit Connection"));
    await session.step(60, "Given Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree")));
    await session.step(61, "When user right-clicks on Platform---Functions---OpenAPI---BDD-Conn-Swagger---Current-Weather-Data-By-City-Name tree node inside browse tree", () => rightClickOn(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger---Current-Weather-Data-By-City-Name tree node inside browse tree")));
    await session.step(62, "And user picks \"Run\" from the open menu", () => pickFromOpenMenu(page, "Run"));
    await session.step(63, "Then \"Current Weather Data By City Name\" dialog should be visible", () => shouldBe(page, el("\"Current Weather Data By City Name\" dialog"), "visible"));
    await session.step(64, "When user enters \"London\" into Q input in \"Current Weather Data By City Name\" dialog", () => enterInto(page, "London", el("Q input in \"Current Weather Data By City Name\" dialog")));
    await session.step(65, "Given user watches the task bar", () => watchTaskBar(page));
    await session.step(66, "When user clicks on OK button in \"Current Weather Data By City Name\" dialog", () => clickOn(page, el("OK button in \"Current Weather Data By City Name\" dialog")));
    await session.step(67, "Then the task bar should have finished \"Running Current Weather Data By City Name\"", () => taskBarFinished(page, "Running Current Weather Data By City Name"));
    await session.step(68, "And the table should have 1 row", () => rowCount(page, 1));
    await session.step(69, "And the table should have a column \"weather/main\"", () => hasColumn(page, "weather/main"));
    await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
