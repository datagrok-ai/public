/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/functions/function-view-export.feature
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
import {clickOn, downloadThrough, fileMatchingDownloaded, pressKey, shouldBe, shouldHaveValue, shouldNotBe, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {scriptOnServer, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The function view exports a run to Excel, PowerPoint or the script's own formats", () => {
  const session = feature(test, "features/functions/function-view-export.feature", import.meta.url);
  test("The function view exports a run to Excel, PowerPoint or the script's own formats", {tag: ["@journey", "@realizes:views.functions"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And the runs of the \"BDDFuncExport{time}\" function are deleted at the end", () => runsDeleted(page, session.text("BDDFuncExport{time}")));
    await session.step(15, "And the runs of the \"BDDFuncOwnExport{time}\" function are deleted at the end", () => runsDeleted(page, session.text("BDDFuncOwnExport{time}")));
    await session.step(16, "And a script \"BDDFuncExport{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDFuncExport{time}"), session.text("//language: javascript\n//input: int rows = 5 {caption: Rows}\n//output: dataframe curve {caption: Curve; viewer: Line chart(x: \"x\", y: \"y\") | Grid}\n//output: int count {caption: Count}\n//meta.customExports: [{\"name\": \"Summary for the lab\", \"function\": \"Admin:BDDFuncOwnExport{time}\"}]\nconst xs = Array.from({length: rows}, (_, i) => i);\ncurve = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'x', xs), DG.Column.fromList('int', 'y', xs.map((x) => x * x))]);\ncount = rows;")));
    await session.step(27, "And a script \"BDDFuncOwnExport{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDFuncOwnExport{time}"), session.text("//language: javascript\n//input: int rows = 5 {caption: Rows}\n//output: int count {caption: Count}\n//meta.features: {\"export\": false}\n//meta.customExports: [{\"name\": \"Summary for the lab\", \"function\": \"Admin:BDDFuncExport{time}\"}]\ncount = rows;")));
    await session.step(36, "And user watches downloads", () => watchDownloads(page));
    await run.scenario("The download icon offers Excel, PowerPoint and then the script's own exports", async () => {
      await session.step(39, "When user opens the address \"/func/Admin.BDDFuncExport{time}\"", () => openAddress(page, session.text("/func/Admin.BDDFuncExport{time}")));
      await session.step(40, "And user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(41, "Then \"Line chart\" viewer should be visible", () => shouldBe(page, el("\"Line chart\" viewer"), "visible"));
      await session.step(42, "When user clicks on \"Export the results\" icon", () => clickOn(page, el("\"Export the results\" icon")));
      await session.step(43, "Then \"Excel\" menu item should be visible", () => shouldBe(page, el("\"Excel\" menu item"), "visible"));
      await session.step(44, "And \"PowerPoint\" menu item should be visible", () => shouldBe(page, el("\"PowerPoint\" menu item"), "visible"));
      await session.step(45, "And \"Summary for the lab\" menu item should be visible", () => shouldBe(page, el("\"Summary for the lab\" menu item"), "visible"));
      await session.step(46, "When user presses Escape", () => pressKey(page, "Escape"));
    });
    await run.scenario("Excel downloads the run as a workbook", async () => {
      await session.step(49, "When user clicks on \"Export the results\" icon", () => clickOn(page, el("\"Export the results\" icon")));
      await session.step(50, "And user downloads a file through \"Excel\" menu item", () => downloadThrough(page, el("\"Excel\" menu item")));
      await session.step(51, "Then a file matching \"\\.xlsx$\" should have been downloaded", () => fileMatchingDownloaded(page, "\\.xlsx$"));
    });
    await run.scenario("PowerPoint downloads the run as a presentation", async () => {
      await session.step(54, "When user clicks on \"Export the results\" icon", () => clickOn(page, el("\"Export the results\" icon")));
      await session.step(55, "And user downloads a file through \"PowerPoint\" menu item", () => downloadThrough(page, el("\"PowerPoint\" menu item")));
      await session.step(56, "Then a file matching \"\\.pptx$\" should have been downloaded", () => fileMatchingDownloaded(page, "\\.pptx$"));
    });
    await run.scenario("A script that turns the built-in exports off offers only its own", async () => {
      await session.step(59, "When user opens the address \"/func/Admin.BDDFuncOwnExport{time}\"", () => openAddress(page, session.text("/func/Admin.BDDFuncOwnExport{time}")));
      await session.step(60, "And user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(61, "Then \"Count\" input should have value \"5\"", () => shouldHaveValue(page, el("\"Count\" input"), "5"));
      await session.step(62, "And \"Export to Summary for the lab\" icon should be visible", () => shouldBe(page, el("\"Export to Summary for the lab\" icon"), "visible"));
      await session.step(63, "And \"Export the results\" icon should not be visible", () => shouldNotBe(page, el("\"Export the results\" icon"), "visible"));
    });
    run.finish();
  });
});
