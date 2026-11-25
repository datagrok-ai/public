/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/flow-dashboards.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, views.space, public#4675, public#4676]
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
import {clickDialogLink} from '../../bindings/close-tabs.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, noProjectOnServer, openApp, openProjectWithTable, saveFromDialog, savedWithDataSync, scriptOnServer, simpleModeOff, spaceOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The dashboards a flow publishes", () => {
  const session = feature(test, "features/apps/flow-dashboards.feature", import.meta.url);
  test("A flow published as a new dashboard leaves the first one with its table", {tag: ["@realizes:views.projects", "@realizes:views.space", "@realizes:public#4675"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "Given no project named \"BDDFlowFirst{time}, BDDFlowSecond{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDFlowFirst{time}, BDDFlowSecond{time}")));
    await session.step(19, "And a script \"BDDFlowPublish{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDFlowPublish{time}"), session.text("//language: flow\n//tags: flow\n//output: dataframe result\n{\"version\":\"2.0\",\"name\":\"BDDFlowPublish{time}\",\"description\":\"\",\"author\":\"\",\"created\":\"\",\"modified\":\"\",\"nodes\":[{\"id\":\"a\",\"typeName\":\"DG Functions/Uncategorized/OpenFile\",\"label\":\"Open File: cars.csv\",\"description\":\"\",\"pos\":{\"x\":40,\"y\":40},\"properties\":{},\"inputValues\":{\"fullPath\":\"System:DemoFiles/cars.csv\",\"sheetName\":\"\"}},{\"id\":\"b\",\"typeName\":\"Outputs/Table Output\",\"label\":\"Table Output\",\"description\":\"\",\"pos\":{\"x\":400,\"y\":40},\"properties\":{\"paramName\":\"result\"},\"inputValues\":{}}],\"connections\":[{\"id\":\"c\",\"source\":\"a\",\"sourceOutput\":\"result\",\"target\":\"b\",\"targetInput\":\"table\"}],\"annotations\":[],\"groups\":[],\"metadata\":{\"settings\":{\"scriptName\":\"BDDFlowPublish{time}\",\"scriptDescription\":\"\",\"tags\":[\"flow\"]}}}")));
    await session.step(26, "And user opens the \"Flow\" app", () => openApp(page, "Flow"));
    await session.step(27, "When user clicks on \"BDDFlowPublish{time}\" recent flow", () => clickOn(page, el(session.text("\"BDDFlowPublish{time}\" recent flow"))));
    await session.step(28, "Then the \"nodes\" reading of flow editor widget should be 2", () => readingIs(page, "nodes", el("flow editor widget"), 2));
    await session.step(29, "When user clicks on flow run button", () => clickOn(page, el("flow run button")));
    await session.step(30, "Then the \"computed tables\" reading of flow editor widget should be 1", () => readingIs(page, "computed tables", el("flow editor widget"), 1));
    await session.step(31, "When user clicks on flow save button", () => clickOn(page, el("flow save button")));
    await session.step(32, "And user clicks on \"Save & publish…\" button in \"Save Flow\" dialog", () => clickOn(page, el("\"Save & publish…\" button in \"Save Flow\" dialog")));
    await session.step(33, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(34, "When user enters \"BDDFlowFirst{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDFlowFirst{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(35, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
    await session.step(36, "And user clicks on flow save button", () => clickOn(page, el("flow save button")));
    await session.step(37, "And user selects \"As a new dashboard\" in Publish choice input in \"Save Flow\" dialog", () => selectIn(page, "As a new dashboard", el("Publish choice input in \"Save Flow\" dialog")));
    await session.step(38, "And user clicks on \"Save & publish…\" button in \"Save Flow\" dialog", () => clickOn(page, el("\"Save & publish…\" button in \"Save Flow\" dialog")));
    await session.step(39, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(40, "When user enters \"BDDFlowSecond{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDFlowSecond{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(41, "And user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
    await session.step(42, "Then the \"cars\" table of the \"BDDFlowSecond{time}\" project should be saved with data sync", () => savedWithDataSync(page, "cars", session.text("BDDFlowSecond{time}")));
    await session.step(43, "And the \"cars\" table of the \"BDDFlowFirst{time}\" project should be saved with data sync", () => savedWithDataSync(page, "cars", session.text("BDDFlowFirst{time}")));
    await session.step(44, "When user closes all views", () => closeAllViews(page));
    await session.step(45, "And user opens the \"BDDFlowFirst{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDDFlowFirst{time}")));
    await session.step(46, "Then the table should have 30 rows", () => rowCount(page, 30));
  });
  test("A dashboard published while its flow is saved into a space opens with its data", {tag: ["@realizes:views.projects", "@realizes:views.space", "@realizes:public#4676"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(50, "Given no project named \"BDDFlowMoved{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDFlowMoved{time}")));
    await session.step(51, "And a script \"BDDFlowMoved{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDFlowMoved{time}"), session.text("//language: flow\n//tags: flow\n//output: dataframe result\n{\"version\":\"2.0\",\"name\":\"BDDFlowMoved{time}\",\"description\":\"\",\"author\":\"\",\"created\":\"\",\"modified\":\"\",\"nodes\":[{\"id\":\"a\",\"typeName\":\"DG Functions/Uncategorized/OpenFile\",\"label\":\"Open File: cars.csv\",\"description\":\"\",\"pos\":{\"x\":40,\"y\":40},\"properties\":{},\"inputValues\":{\"fullPath\":\"System:DemoFiles/cars.csv\",\"sheetName\":\"\"}},{\"id\":\"b\",\"typeName\":\"Outputs/Table Output\",\"label\":\"Table Output\",\"description\":\"\",\"pos\":{\"x\":400,\"y\":40},\"properties\":{\"paramName\":\"result\"},\"inputValues\":{}}],\"connections\":[{\"id\":\"c\",\"source\":\"a\",\"sourceOutput\":\"result\",\"target\":\"b\",\"targetInput\":\"table\"}],\"annotations\":[],\"groups\":[],\"metadata\":{\"settings\":{\"scriptName\":\"BDDFlowMoved{time}\",\"scriptDescription\":\"\",\"tags\":[\"flow\"]}}}")));
    await session.step(58, "And a space named \"BDDFlowHome{time}\" is on the server", () => spaceOnServer(page, session.text("BDDFlowHome{time}")));
    await session.step(59, "And user opens the \"Flow\" app", () => openApp(page, "Flow"));
    await session.step(60, "When user clicks on \"BDDFlowMoved{time}\" recent flow", () => clickOn(page, el(session.text("\"BDDFlowMoved{time}\" recent flow"))));
    await session.step(61, "Then the \"nodes\" reading of flow editor widget should be 2", () => readingIs(page, "nodes", el("flow editor widget"), 2));
    await session.step(62, "When user clicks on flow run button", () => clickOn(page, el("flow run button")));
    await session.step(63, "Then the \"computed tables\" reading of flow editor widget should be 1", () => readingIs(page, "computed tables", el("flow editor widget"), 1));
    await session.step(64, "When user clicks on flow save button", () => clickOn(page, el("flow save button")));
    await session.step(65, "And user clicks on the \"Choose…\" link in \"Save Flow\" dialog", () => clickDialogLink(page, "Choose…", "Save Flow"));
    await session.step(66, "And user clicks on \"BDDFlowHome{time}\" tree node in \"Save Flow\" dialog", () => clickOn(page, el(session.text("\"BDDFlowHome{time}\" tree node in \"Save Flow\" dialog"))));
    await session.step(67, "And user clicks on \"Save & publish…\" button in \"Save Flow\" dialog", () => clickOn(page, el("\"Save & publish…\" button in \"Save Flow\" dialog")));
    await session.step(68, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(69, "And Name text input in \"Save project\" dialog should have value \"BDDFlowMoved{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("BDDFlowMoved{time}")));
    await session.step(70, "When user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
    await session.step(71, "And user closes all views", () => closeAllViews(page));
    await session.step(72, "And user opens the \"BDDFlowMoved{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDDFlowMoved{time}")));
    await session.step(73, "Then \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
    await session.step(74, "And the table should have 30 rows", () => rowCount(page, 30));
  });
  test("So does one published into the same space under the flow's name", {tag: ["@realizes:views.projects", "@realizes:views.space", "@realizes:public#4676"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(78, "Given no project named \"BDDFlowInSpace{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDFlowInSpace{time}")));
    await session.step(79, "And a script \"BDDFlowInSpace{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDFlowInSpace{time}"), session.text("//language: flow\n//tags: flow\n//output: dataframe result\n{\"version\":\"2.0\",\"name\":\"BDDFlowInSpace{time}\",\"description\":\"\",\"author\":\"\",\"created\":\"\",\"modified\":\"\",\"nodes\":[{\"id\":\"a\",\"typeName\":\"DG Functions/Uncategorized/OpenFile\",\"label\":\"Open File: cars.csv\",\"description\":\"\",\"pos\":{\"x\":40,\"y\":40},\"properties\":{},\"inputValues\":{\"fullPath\":\"System:DemoFiles/cars.csv\",\"sheetName\":\"\"}},{\"id\":\"b\",\"typeName\":\"Outputs/Table Output\",\"label\":\"Table Output\",\"description\":\"\",\"pos\":{\"x\":400,\"y\":40},\"properties\":{\"paramName\":\"result\"},\"inputValues\":{}}],\"connections\":[{\"id\":\"c\",\"source\":\"a\",\"sourceOutput\":\"result\",\"target\":\"b\",\"targetInput\":\"table\"}],\"annotations\":[],\"groups\":[],\"metadata\":{\"settings\":{\"scriptName\":\"BDDFlowInSpace{time}\",\"scriptDescription\":\"\",\"tags\":[\"flow\"]}}}")));
    await session.step(86, "And a space named \"BDDFlowSpace{time}\" is on the server", () => spaceOnServer(page, session.text("BDDFlowSpace{time}")));
    await session.step(87, "And user opens the \"Flow\" app", () => openApp(page, "Flow"));
    await session.step(88, "When user clicks on \"BDDFlowInSpace{time}\" recent flow", () => clickOn(page, el(session.text("\"BDDFlowInSpace{time}\" recent flow"))));
    await session.step(89, "Then the \"nodes\" reading of flow editor widget should be 2", () => readingIs(page, "nodes", el("flow editor widget"), 2));
    await session.step(90, "When user clicks on flow run button", () => clickOn(page, el("flow run button")));
    await session.step(91, "Then the \"computed tables\" reading of flow editor widget should be 1", () => readingIs(page, "computed tables", el("flow editor widget"), 1));
    await session.step(92, "When user clicks on flow save button", () => clickOn(page, el("flow save button")));
    await session.step(93, "And user clicks on the \"Choose…\" link in \"Save Flow\" dialog", () => clickDialogLink(page, "Choose…", "Save Flow"));
    await session.step(94, "And user clicks on \"BDDFlowSpace{time}\" tree node in \"Save Flow\" dialog", () => clickOn(page, el(session.text("\"BDDFlowSpace{time}\" tree node in \"Save Flow\" dialog"))));
    await session.step(95, "And user clicks on \"Save & publish…\" button in \"Save Flow\" dialog", () => clickOn(page, el("\"Save & publish…\" button in \"Save Flow\" dialog")));
    await session.step(96, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(97, "And Name text input in \"Save project\" dialog should have value \"BDDFlowInSpace{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("BDDFlowInSpace{time}")));
    await session.step(98, "When user clicks on space input in \"Save project\" dialog", () => clickOn(page, el("space input in \"Save project\" dialog")));
    await session.step(99, "And user clicks on \"BDDFlowSpace{time}\" tree node", () => clickOn(page, el(session.text("\"BDDFlowSpace{time}\" tree node"))));
    await session.step(100, "Then space input in \"Save project\" dialog should contain text \"BDDFlowSpace{time}\"", () => shouldContainText(page, el("space input in \"Save project\" dialog"), session.text("BDDFlowSpace{time}")));
    await session.step(101, "When user clicks on OK in the Save project dialog and the project uploads", () => saveFromDialog(page));
    await session.step(102, "And user closes all views", () => closeAllViews(page));
    await session.step(103, "And user opens the \"BDDFlowInSpace{time}\" project and waits for its table", () => openProjectWithTable(page, session.text("BDDFlowInSpace{time}")));
    await session.step(104, "Then \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
    await session.step(105, "And the table should have 30 rows", () => rowCount(page, 30));
  });
});
