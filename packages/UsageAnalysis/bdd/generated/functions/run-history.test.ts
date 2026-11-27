/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/functions/run-history.feature
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
import {ownRunsDeleted, runsDeleted, spaceLinksRun} from '../../bindings/function-view.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, hoverOver, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, typeInto, visibleCount} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {dialogCloses, openDataset, packageInstalled, scriptOnServer, simpleModeOff, spaceOnServer, urlShouldContain} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Run history records every run of the function view and brings any back", () => {
  const session = feature(test, "features/functions/run-history.feature", import.meta.url);
  test("Run history records every run of the function view and brings any back", {tag: ["@journey", "@realizes:views.functions"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And the runs of the \"BDDRunHistory{time}\" function are deleted at the end", () => runsDeleted(page, session.text("BDDRunHistory{time}")));
    await session.step(19, "And the runs of the \"BDDRunLive{time}\" function are deleted at the end", () => runsDeleted(page, session.text("BDDRunLive{time}")));
    await session.step(20, "And a script \"BDDRunHistory{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDRunHistory{time}"), "//language: javascript\n//input: int rows = 5 {caption: Rows}\n//output: int total {caption: Total}\ntotal = rows * 10;"));
    await session.step(27, "And a script \"BDDRunLive{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDRunLive{time}"), "//language: javascript\n//input: int rows = 5 {caption: Rows}\n//output: int total {caption: Total}\n//meta.runOnInput: true\ntotal = rows * 10;"));
    await session.step(35, "And a space named \"BDDRunSpace{time}\" is on the server", () => spaceOnServer(page, session.text("BDDRunSpace{time}")));
    await run.scenario("Every RUN is recorded with its own number", async () => {
      await session.step(38, "When user opens the address \"/func/Admin.BDDRunHistory{time}?rows=5\"", () => openAddress(page, session.text("/func/Admin.BDDRunHistory{time}?rows=5")));
      await session.step(39, "And user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(40, "Then run header should contain text \"Run #1\"", () => shouldContainText(page, el("run header"), "Run #1"));
      await session.step(41, "When user enters \"7\" into \"Rows\" input", () => enterInto(page, "7", el("\"Rows\" input")));
      await session.step(42, "And user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(43, "Then run header should contain text \"Run #2\"", () => shouldContainText(page, el("run header"), "Run #2"));
      await session.step(44, "And \"Total\" input should have value \"70\"", () => shouldHaveValue(page, el("\"Total\" input"), "70"));
    });
    await run.scenario("The history icon docks the runs, each with what changed against the run before it", async () => {
      await session.step(47, "When user clicks on \"Run history\" icon", () => clickOn(page, el("\"Run history\" icon")));
      await session.step(48, "Then run history should be visible", () => shouldBe(page, el("run history"), "visible"));
      await session.step(49, "And \"Run #2\" run should be visible", () => shouldBe(page, el("\"Run #2\" run"), "visible"));
      await session.step(50, "And \"Run #1\" run should be visible", () => shouldBe(page, el("\"Run #1\" run"), "visible"));
      await session.step(51, "And changes of \"Run #2\" run should contain text \"rows 5 → 7\"", () => shouldContainText(page, el("changes of \"Run #2\" run"), "rows 5 → 7"));
      await session.step(52, "And changes of \"Run #2\" run should contain text \"total = 70\"", () => shouldContainText(page, el("changes of \"Run #2\" run"), "total = 70"));
    });
    await run.scenario("The pencil names a run, and the name shows in the header and in the history", async () => {
      await session.step(55, "When user clicks on \"Edit title, description and tags\" icon in run header", () => clickOn(page, el("\"Edit title, description and tags\" icon in run header")));
      await session.step(56, "Then \"Run #2\" dialog should be visible", () => shouldBe(page, el("\"Run #2\" dialog"), "visible"));
      await session.step(57, "When user enters \"Baseline\" into \"Title\" input in \"Run #2\" dialog", () => enterInto(page, "Baseline", el("\"Title\" input in \"Run #2\" dialog")));
      await session.step(58, "And user clicks on OK button in \"Run #2\" dialog", () => clickOn(page, el("OK button in \"Run #2\" dialog")));
      await session.step(59, "Then run header should contain text \"Baseline\"", () => shouldContainText(page, el("run header"), "Baseline"));
      await session.step(60, "And \"Baseline\" run should be visible", () => shouldBe(page, el("\"Baseline\" run"), "visible"));
    });
    await run.scenario("The star keeps a run, and the Starred filter shows only the starred ones", async () => {
      await session.step(63, "When user clicks on run header star", () => clickOn(page, el("run header star")));
      await session.step(64, "And user clicks on \"Starred\" run filter", () => clickOn(page, el("\"Starred\" run filter")));
      await session.step(65, "Then there should be 1 visible run in run history", () => visibleCount(page, 1, el("run in run history")));
      await session.step(66, "And \"Baseline\" run should be visible", () => shouldBe(page, el("\"Baseline\" run"), "visible"));
      await session.step(67, "When user clicks on \"Starred\" run filter", () => clickOn(page, el("\"Starred\" run filter")));
      await session.step(68, "Then there should be 2 visible run in run history", () => visibleCount(page, 2, el("run in run history")));
    });
    await run.scenario("The search finds a run by its name", async () => {
      await session.step(71, "When user types \"Baseline\" into run history search", () => typeInto(page, "Baseline", el("run history search")));
      await session.step(72, "Then there should be 1 visible run in run history", () => visibleCount(page, 1, el("run in run history")));
      await session.step(73, "When user clears run history search", () => clearField(page, el("run history search")));
      await session.step(74, "Then there should be 2 visible run in run history", () => visibleCount(page, 2, el("run in run history")));
    });
    await run.scenario("A click on a run brings its inputs back at the run's own address", async () => {
      await session.step(77, "When user clicks on \"Run #1\" run", () => clickOn(page, el("\"Run #1\" run")));
      await session.step(78, "Then \"Rows\" input should have value \"5\"", () => shouldHaveValue(page, el("\"Rows\" input"), "5"));
      await session.step(79, "And run header should contain text \"Run #1\"", () => shouldContainText(page, el("run header"), "Run #1"));
      await session.step(80, "And the page address should contain \"/runs/1\"", () => urlShouldContain(page, "/runs/1"));
    });
    await run.scenario("The author deletes a run of their own", async () => {
      await session.step(83, "When user hovers over \"Run #1\" run", () => hoverOver(page, el("\"Run #1\" run")));
      await session.step(84, "And user clicks on delete of \"Run #1\" run", () => clickOn(page, el("delete of \"Run #1\" run")));
      await session.step(85, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(86, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(87, "Then \"Run #1\" run should not be visible", () => shouldNotBe(page, el("\"Run #1\" run"), "visible"));
      await session.step(88, "And \"Baseline\" run should be visible", () => shouldBe(page, el("\"Baseline\" run"), "visible"));
    });
    await run.scenario("The panel leaves with its view and comes back with it", async () => {
      await session.step(91, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(92, "Then run history should not be visible", () => shouldNotBe(page, el("run history"), "visible"));
      await session.step(93, "When user clicks on \"BDDRunHistory{time}\" view tab", () => clickOn(page, el(session.text("\"BDDRunHistory{time}\" view tab"))));
      await session.step(94, "Then run history should be visible", () => shouldBe(page, el("run history"), "visible"));
    });
    await run.scenario("Save to space links the run into the space picked", async () => {
      await session.step(97, "Given the \"PowerPack\" package is installed", () => packageInstalled(page, "PowerPack"));
      await session.step(98, "When user clicks on \"BDDRunHistory{time}\" view tab", () => clickOn(page, el(session.text("\"BDDRunHistory{time}\" view tab"))));
      await session.step(99, "And user clicks on \"Baseline\" run", () => clickOn(page, el("\"Baseline\" run")));
      await session.step(100, "Then run header should contain text \"Baseline\"", () => shouldContainText(page, el("run header"), "Baseline"));
      await session.step(101, "When user clicks on \"Save this run to a space\" icon in run header", () => clickOn(page, el("\"Save this run to a space\" icon in run header")));
      await session.step(102, "Then \"Save Baseline to a space\" dialog should be visible", () => shouldBe(page, el("\"Save Baseline to a space\" dialog"), "visible"));
      await session.step(103, "And space input in \"Save Baseline to a space\" dialog should contain text \"Your own projects\"", () => shouldContainText(page, el("space input in \"Save Baseline to a space\" dialog"), "Your own projects"));
      await session.step(104, "When user clicks on space input in \"Save Baseline to a space\" dialog", () => clickOn(page, el("space input in \"Save Baseline to a space\" dialog")));
      await session.step(105, "And user clicks on \"BDDRunSpace{time}\" tree node", () => clickOn(page, el(session.text("\"BDDRunSpace{time}\" tree node"))));
      await session.step(106, "And user clicks on OK button in \"Save Baseline to a space\" dialog", () => clickOn(page, el("OK button in \"Save Baseline to a space\" dialog")));
      await session.step(107, "Then the \"Save Baseline to a space\" dialog should close", () => dialogCloses(page, "Save Baseline to a space"));
      await session.step(108, "And the space \"BDDRunSpace{time}\" on the server should link a run of the \"BDDRunHistory{time}\" function", () => spaceLinksRun(page, session.text("BDDRunSpace{time}"), session.text("BDDRunHistory{time}")));
    });
    await run.scenario("A function that reruns on input records one run once the inputs rest", async () => {
      await session.step(111, "When user opens the address \"/func/Admin.BDDRunLive{time}?rows=5\"", () => openAddress(page, session.text("/func/Admin.BDDRunLive{time}?rows=5")));
      await session.step(112, "And user enters \"6\" into \"Rows\" input", () => enterInto(page, "6", el("\"Rows\" input")));
      await session.step(113, "And user enters \"8\" into \"Rows\" input", () => enterInto(page, "8", el("\"Rows\" input")));
      await session.step(114, "Then \"Total\" input should have value \"80\"", () => shouldHaveValue(page, el("\"Total\" input"), "80"));
      await session.step(115, "And run header should contain text \"Run #\"", () => shouldContainText(page, el("run header"), "Run #"));
      await session.step(116, "When user clicks on \"Run history\" icon", () => clickOn(page, el("\"Run history\" icon")));
      await session.step(117, "Then \"Run #1\" run should be visible", () => shouldBe(page, el("\"Run #1\" run"), "visible"));
      await session.step(118, "And changes of \"Run #1\" run should contain text \"rows 8\"", () => shouldContainText(page, el("changes of \"Run #1\" run"), "rows 8"));
    });
    await run.scenario("A package function's history lists its runs", async () => {
      await session.step(121, "Given the \"Compute2\" package is installed", () => packageInstalled(page, "Compute2"));
      await session.step(122, "And the runs of the \"ObjectCooling2\" function made from now on are deleted at the end", () => ownRunsDeleted(page, "ObjectCooling2"));
      await session.step(123, "When user opens the address \"/func/Compute2.ObjectCooling2\"", () => openAddress(page, "/func/Compute2.ObjectCooling2"));
      await session.step(124, "And user clicks on RUN button", () => clickOn(page, el("RUN button")));
      await session.step(125, "Then run header should contain text \"Run #\"", () => shouldContainText(page, el("run header"), "Run #"));
      await session.step(126, "When user clicks on \"Run history\" icon", () => clickOn(page, el("\"Run history\" icon")));
      await session.step(127, "Then run history should be visible", () => shouldBe(page, el("run history"), "visible"));
      await session.step(128, "And first run in run history should be visible", () => shouldBe(page, el("first run in run history"), "visible"));
    });
    run.finish();
  });
});
