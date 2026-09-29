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
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, followingShouldBe, isExpanded, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, connectionsOnServer, noConnectionOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, noErrors, openContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Importing an OpenAPI (Swagger) file as a connection", () => {
  const session = feature(test, "features/connections/connections-import-swagger.feature", import.meta.url);
  test("The imported file is a connection with a query per operation", {tag: ["@connections", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And no connection named \"BDD-Conn-Swagger\" is on the server", () => noConnectionOnServer(page, "BDD-Conn-Swagger"));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "When user uploads \"fixtures/bdd-swagger.yaml\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/bdd-swagger.yaml", el("\"Open local file\" icon inside browse toolbar")));
    await session.step(25, "Then 1 connection named \"BDD-Conn-Swagger\" should be on the server", () => connectionsOnServer(page, 1, "BDD-Conn-Swagger"));
    await session.step(26, "And the \"bdd-swagger.yaml\" view should be current", () => viewIsCurrent(page, "bdd-swagger.yaml"));
    await session.step(27, "When user closes all views", () => closeAllViews(page));
    await session.step(29, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(30, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(31, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(32, "And Platform---Functions---OpenAPI tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions---OpenAPI tree node inside browse tree")));
    await session.step(35, "Then Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree should be visible", () => shouldBe(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree"), "visible"));
    await session.step(36, "Given Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree")));
    await session.step(37, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Current-Weather-Data-By-City-Name tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Cities-Within-a-Rectangle-Zone tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Cities-In-Cycle tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---5-day/3-hour-Forecast-By-City-Name tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Call-Current-UV-Data-For-One-Location-By-Geographic-Coordinates tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Call-Forecast-UV-Data-For-One-Location-By-Geographic-Coordinates tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Call-Historical-UV-Data-For-One-Location tree node inside browse tree"]]), [["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Current-Weather-Data-By-City-Name tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Cities-Within-a-Rectangle-Zone tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Cities-In-Cycle tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---5-day/3-hour-Forecast-By-City-Name tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Call-Current-UV-Data-For-One-Location-By-Geographic-Coordinates tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Call-Forecast-UV-Data-For-One-Location-By-Geographic-Coordinates tree node inside browse tree"],["Platform---Functions---OpenAPI---BDD-Conn-Swagger---Call-Historical-UV-Data-For-One-Location tree node inside browse tree"]]);
    await session.step(45, "When user opens the context menu of Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree", () => openContextMenu(page, el("Platform---Functions---OpenAPI---BDD-Conn-Swagger tree node inside browse tree")));
    await session.step(46, "Then the open menu should list \"Edit...\"", () => menuLists(page, "Edit..."));
    await session.step(47, "And the open menu should list \"Test connection\"", () => menuLists(page, "Test connection"));
    await session.step(48, "When user closes the context menu", () => closeContextMenu(page));
    await session.step(49, "Then no errors should have been logged", () => noErrors(page));
  });
});
