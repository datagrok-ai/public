/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/functions/function-view.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.functions]
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
import {runsDeleted} from '../../bindings/function-view.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, clipboardContains, enterInto, hoverOver, pressKey, pressKeyIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {dialogCloses, scriptOnServer, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The function view shows a run's results and keeps them in step with the form", () => {
  const session = feature(test, "features/functions/function-view.feature", import.meta.url);
  test("The function view shows a run's results and keeps them in step with the form", {tag: ["@journey", "@realizes:views.functions"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And the runs of the \"BDDFuncView{time}\" function are deleted at the end", () => runsDeleted(page, session.text("BDDFuncView{time}")));
    await session.step(17, "And a script \"BDDFuncView{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDFuncView{time}"), "//language: javascript\n//input: int rows = 5 {caption: Rows; nullable: false}\n//input: double slope = 2 {caption: Slope}\n//output: dataframe curve {caption: Curve; viewer: Line chart(x: \"x\", y: \"y\", title: \"Response\") | Grid}\n//output: double total {caption: Total}\nconst xs = Array.from({length: rows}, (_, i) => i);\ncurve = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'x', xs), DG.Column.fromList('double', 'y', xs.map((x) => slope * x))]);\ntotal = xs.reduce((s, x) => s + slope * x, 0);"));
    await run.scenario("RUN shows the results with the viewers the script declares", async () => {
      await session.step(30, "When user opens the address \"/func/Admin.BDDFuncView{time}\"", () => openAddress(page, session.text("/func/Admin.BDDFuncView{time}")));
      await session.step(31, "Then \"Rows\" input should have value \"5\"", () => shouldHaveValue(page, el("\"Rows\" input"), "5"));
      await session.step(32, "When user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(33, "Then run header should contain text \"Run #1\"", () => shouldContainText(page, el("run header"), "Run #1"));
      await session.step(34, "And \"Line chart\" viewer should be visible", () => shouldBe(page, el("\"Line chart\" viewer"), "visible"));
      await session.step(35, "And function results should contain text \"Response\"", () => shouldContainText(page, el("function results"), "Response"));
      await session.step(36, "And \"Grid\" viewer should be visible", () => shouldBe(page, el("\"Grid\" viewer"), "visible"));
      await session.step(37, "And \"Total\" input should have value \"20\"", () => shouldHaveValue(page, el("\"Total\" input"), "20"));
    });
    await run.scenario("A changed input marks the results as outdated until the next RUN", async () => {
      await session.step(40, "When user enters \"7\" into \"Rows\" input", () => enterInto(page, "7", el("\"Rows\" input")));
      await session.step(41, "Then function results should contain text \"Inputs changed\"", () => shouldContainText(page, el("function results"), "Inputs changed"));
      await session.step(42, "When user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(43, "Then function results should not contain text \"Inputs changed\"", () => shouldNotContainText(page, el("function results"), "Inputs changed"));
      await session.step(44, "And \"Total\" input should have value \"42\"", () => shouldHaveValue(page, el("\"Total\" input"), "42"));
      await session.step(45, "And run header should contain text \"Run #2\"", () => shouldContainText(page, el("run header"), "Run #2"));
    });
    await run.scenario("A re-run keeps a viewer the user closed closed, and Ctrl+Z brings it back", async () => {
      await session.step(48, "When user hovers over line chart viewer", () => hoverOver(page, el("line chart viewer")));
      await session.step(49, "And user clicks on close icon of line chart viewer", () => clickOn(page, el("close icon of line chart viewer")));
      await session.step(50, "Then \"Line chart\" viewer should not be visible", () => shouldNotBe(page, el("\"Line chart\" viewer"), "visible"));
      await session.step(51, "When user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(52, "Then run header should contain text \"Run #3\"", () => shouldContainText(page, el("run header"), "Run #3"));
      await session.step(53, "And \"Line chart\" viewer should not be visible", () => shouldNotBe(page, el("\"Line chart\" viewer"), "visible"));
      await session.step(54, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(55, "Then \"Line chart\" viewer should be visible", () => shouldBe(page, el("\"Line chart\" viewer"), "visible"));
    });
    await run.scenario("Enter runs only a form whose inputs are valid", async () => {
      await session.step(58, "When user clears \"Rows\" input", () => clearField(page, el("\"Rows\" input")));
      await session.step(59, "And user presses Enter in \"Rows\" input", () => pressKeyIn(page, "Enter", el("\"Rows\" input")));
      await session.step(60, "Then run header should contain text \"Run #3\"", () => shouldContainText(page, el("run header"), "Run #3"));
      await session.step(61, "When user enters \"3\" into \"Rows\" input", () => enterInto(page, "3", el("\"Rows\" input")));
      await session.step(62, "And user presses Enter in \"Rows\" input", () => pressKeyIn(page, "Enter", el("\"Rows\" input")));
      await session.step(63, "Then run header should contain text \"Run #4\"", () => shouldContainText(page, el("run header"), "Run #4"));
      await session.step(64, "And \"Total\" input should have value \"6\"", () => shouldHaveValue(page, el("\"Total\" input"), "6"));
    });
    await run.scenario("The copy icon gives a link with the current inputs, and a link opens the form with its values", async () => {
      await session.step(67, "When user clicks on copy icon", () => clickOn(page, el("copy icon")));
      await session.step(68, "Then the clipboard should contain text \"/func/Admin.BDDFuncView{time}?rows=3\"", () => clipboardContains(page, session.text("/func/Admin.BDDFuncView{time}?rows=3")));
      await session.step(69, "When user opens the address \"/func/Admin.BDDFuncView{time}?rows=9&slope=1\"", () => openAddress(page, session.text("/func/Admin.BDDFuncView{time}?rows=9&slope=1")));
      await session.step(70, "Then \"Rows\" input should have value \"9\"", () => shouldHaveValue(page, el("\"Rows\" input"), "9"));
      await session.step(71, "And \"Slope\" input should have value \"1\"", () => shouldHaveValue(page, el("\"Slope\" input"), "1"));
    });
    await run.scenario("A built-in function opens by its bare name", async () => {
      await session.step(74, "When user opens the address \"/func/Sin?x=1\"", () => openAddress(page, "/func/Sin?x=1"));
      await session.step(75, "Then \"X\" input should have value \"1\"", () => shouldHaveValue(page, el("\"X\" input"), "1"));
      await session.step(76, "And RUN button should be visible", () => shouldBe(page, el("RUN button"), "visible"));
    });
    await run.scenario("The ribbon adds the result table to the workspace and offers to save the run as a project", async () => {
      await session.step(79, "When user opens the address \"/func/Admin.BDDFuncView{time}?rows=5\"", () => openAddress(page, session.text("/func/Admin.BDDFuncView{time}?rows=5")));
      await session.step(80, "And user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(81, "Then run header should contain text \"Run #\"", () => shouldContainText(page, el("run header"), "Run #"));
      await session.step(82, "When user clicks on \"Save\" button", () => clickOn(page, el("\"Save\" button")));
      await session.step(83, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(84, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(85, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(86, "When user clicks on \"Add table to workspace\" icon", () => clickOn(page, el("\"Add table to workspace\" icon")));
      await session.step(87, "Then the \"curve\" table view should open with 5 rows", () => tableViewOpened(page, "curve", 5));
    });
    run.finish();
  });
});
