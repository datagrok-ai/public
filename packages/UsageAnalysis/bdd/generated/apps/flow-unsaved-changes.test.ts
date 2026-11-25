/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/flow-unsaved-changes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [public#4677]
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
import {clickOn, close, doubleClickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {dialogCloses, openApp, scriptContains, scriptOnServer, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A changed flow asks before it is closed or replaced", () => {
  const session = feature(test, "features/apps/flow-unsaved-changes.feature", import.meta.url);
  test("A changed flow asks before it is closed or replaced", {tag: ["@journey", "@realizes:public#4677"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And a script \"BDDFlowUnsaved{time}\" is on the server:", () => scriptOnServer(page, session.text("BDDFlowUnsaved{time}"), session.text("//language: flow\n//tags: flow\n{\"version\":\"2.0\",\"name\":\"BDDFlowUnsaved{time}\",\"description\":\"\",\"author\":\"\",\"created\":\"\",\"modified\":\"\",\"nodes\":[{\"id\":\"a\",\"typeName\":\"Constants/String\",\"label\":\"String\",\"description\":\"\",\"pos\":{\"x\":100,\"y\":100},\"properties\":{\"value\":\"x\"},\"inputValues\":{}}],\"connections\":[],\"annotations\":[],\"groups\":[],\"metadata\":{\"settings\":{\"scriptName\":\"BDDFlowUnsaved{time}\",\"scriptDescription\":\"\",\"tags\":[\"flow\"]}}}")));
    await session.step(20, "And user opens the \"Flow\" app", () => openApp(page, "Flow"));
    await run.scenario("Closing a changed flow asks, closing an unchanged one does not", async () => {
      await session.step(23, "When user clicks on \"BDDFlowUnsaved{time}\" recent flow", () => clickOn(page, el(session.text("\"BDDFlowUnsaved{time}\" recent flow"))));
      await session.step(24, "Then the \"nodes\" reading of flow editor widget should be 1", () => readingIs(page, "nodes", el("flow editor widget"), 1));
      await session.step(25, "And the \"unsaved changes\" reading of flow editor widget should be \"false\"", () => readingReads(page, "unsaved changes", el("flow editor widget"), "false"));
      await session.step(26, "When user closes \"BDDFlowUnsaved{time}\" view tab", () => close(page, el(session.text("\"BDDFlowUnsaved{time}\" view tab"))));
      await session.step(27, "Then \"BDDFlowUnsaved{time}\" view tab should be absent", () => shouldBe(page, el(session.text("\"BDDFlowUnsaved{time}\" view tab")), "absent"));
      await session.step(28, "Given user opens the \"Flow\" app", () => openApp(page, "Flow"));
      await session.step(29, "When user clicks on \"BDDFlowUnsaved{time}\" recent flow", () => clickOn(page, el(session.text("\"BDDFlowUnsaved{time}\" recent flow"))));
      await session.step(30, "And user clicks on \"Inputs\" pane header in toolbox", () => clickOn(page, el("\"Inputs\" pane header in toolbox")));
      await session.step(31, "And user double-clicks on flow sketcher input item", () => doubleClickOn(page, el("flow sketcher input item")));
      await session.step(32, "Then the \"nodes\" reading of flow editor widget should be 2", () => readingIs(page, "nodes", el("flow editor widget"), 2));
      await session.step(33, "When user closes \"BDDFlowUnsaved{time}\" view tab", () => close(page, el(session.text("\"BDDFlowUnsaved{time}\" view tab"))));
      await session.step(34, "Then \"Close BDDFlowUnsaved{time}?\" dialog should contain text \"Edited the flow.\"", () => shouldContainText(page, el(session.text("\"Close BDDFlowUnsaved{time}?\" dialog")), "Edited the flow."));
      await session.step(35, "When user clicks on CANCEL button in \"Close BDDFlowUnsaved{time}?\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Close BDDFlowUnsaved{time}?\" dialog"))));
      await session.step(36, "Then the \"Close BDDFlowUnsaved{time}?\" dialog should close", () => dialogCloses(page, session.text("Close BDDFlowUnsaved{time}?")));
      await session.step(37, "And the \"nodes\" reading of flow editor widget should be 2", () => readingIs(page, "nodes", el("flow editor widget"), 2));
    });
    await run.scenario("A new flow and another flow ask before they replace it", async () => {
      await session.step(40, "When user picks \"Flow | New…\" from the top menu", () => pickFromTopMenu(page, "Flow | New…"));
      await session.step(41, "Then \"Unsaved changes\" dialog should contain text \"BDDFlowUnsaved{time}\"", () => shouldContainText(page, el("\"Unsaved changes\" dialog"), session.text("BDDFlowUnsaved{time}")));
      await session.step(42, "When user clicks on DISCARD button in \"Unsaved changes\" dialog", () => clickOn(page, el("DISCARD button in \"Unsaved changes\" dialog")));
      await session.step(43, "Then the \"Unsaved changes\" dialog should close", () => dialogCloses(page, "Unsaved changes"));
      await session.step(44, "And the \"nodes\" reading of flow editor widget should be 0", () => readingIs(page, "nodes", el("flow editor widget"), 0));
      await session.step(45, "When user clicks on \"BDDFlowUnsaved{time}\" recent flow", () => clickOn(page, el(session.text("\"BDDFlowUnsaved{time}\" recent flow"))));
      await session.step(46, "Then the \"nodes\" reading of flow editor widget should be 1", () => readingIs(page, "nodes", el("flow editor widget"), 1));
      await session.step(47, "When user double-clicks on flow sketcher input item", () => doubleClickOn(page, el("flow sketcher input item")));
      await session.step(48, "Then the \"nodes\" reading of flow editor widget should be 2", () => readingIs(page, "nodes", el("flow editor widget"), 2));
      await session.step(49, "When user picks \"Flow | Open from platform…\" from the top menu", () => pickFromTopMenu(page, "Flow | Open from platform…"));
      await session.step(50, "And user selects \"BDDFlowUnsaved{time}\" in Flow choice input in \"Open Flow\" dialog", () => selectIn(page, session.text("BDDFlowUnsaved{time}"), el("Flow choice input in \"Open Flow\" dialog")));
      await session.step(51, "And user clicks on OK button in \"Open Flow\" dialog", () => clickOn(page, el("OK button in \"Open Flow\" dialog")));
      await session.step(52, "Then \"Unsaved changes\" dialog should contain text \"BDDFlowUnsaved{time}\"", () => shouldContainText(page, el("\"Unsaved changes\" dialog"), session.text("BDDFlowUnsaved{time}")));
      await session.step(53, "When user clicks on DISCARD button in \"Unsaved changes\" dialog", () => clickOn(page, el("DISCARD button in \"Unsaved changes\" dialog")));
      await session.step(54, "Then the \"nodes\" reading of flow editor widget should be 1", () => readingIs(page, "nodes", el("flow editor widget"), 1));
      await session.step(55, "And the \"unsaved changes\" reading of flow editor widget should be \"false\"", () => readingReads(page, "unsaved changes", el("flow editor widget"), "false"));
    });
    await run.scenario("The close dialog saves the flow first", async () => {
      await session.step(58, "When user double-clicks on flow sketcher input item", () => doubleClickOn(page, el("flow sketcher input item")));
      await session.step(59, "Then the \"nodes\" reading of flow editor widget should be 2", () => readingIs(page, "nodes", el("flow editor widget"), 2));
      await session.step(60, "When user closes \"BDDFlowUnsaved{time}\" view tab", () => close(page, el(session.text("\"BDDFlowUnsaved{time}\" view tab"))));
      await session.step(61, "And user clicks on the \"save the flow\" link in \"Close BDDFlowUnsaved{time}?\" dialog", () => clickDialogLink(page, "save the flow", session.text("Close BDDFlowUnsaved{time}?")));
      await session.step(62, "Then \"BDDFlowUnsaved{time}\" view tab should be absent", () => shouldBe(page, el(session.text("\"BDDFlowUnsaved{time}\" view tab")), "absent"));
      await session.step(63, "And the script \"BDDFlowUnsaved{time}\" on the server should contain \"Inputs/Sketcher Input\"", () => scriptContains(page, session.text("BDDFlowUnsaved{time}"), "Inputs/Sketcher Input"));
    });
    run.finish();
  });
});
