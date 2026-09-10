/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/markup/markup-chrome.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.markup, entities.viewer.action.close-viewer]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKey, shouldBe, shouldHaveText, shouldHaveValue, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, closeContextMenu, hasArea, hasNoArea, menuLists, noErrors, openContextMenu, pickFromContextMenu, propertyShouldBe, readingReads, reportsNoError, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {openViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Markup chrome — Edit content, the title bar and the strip under the content", () => {
  const session = feature(test, "features/viewers/markup/markup-chrome.feature", import.meta.url);
  test("Markup chrome — Edit content, the title bar and the strip under the content", {tag: ["@journey", "@viewers", "@realizes:viewers.markup", "@realizes:entities.viewer.action.close-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "And user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","editable probe"]]));
    await session.step(31, "Then the \"text\" reading of markup viewer should be \"editable probe\"", () => readingReads(page, "text", el("markup viewer"), "editable probe"));
    await session.step(32, "And markup viewer should report no error", () => reportsNoError(page, el("markup viewer")));
    await run.scenario("Edit content... is on the viewer's menu and opens on what is on screen", async () => {
      await session.step(35, "When user opens the context menu of markup viewer", () => openContextMenu(page, el("markup viewer")));
      await session.step(36, "Then the open menu should list \"Edit content...\"", () => menuLists(page, "Edit content..."));
      await session.step(37, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(38, "And user picks \"Edit content...\" from the context menu of markup viewer", () => pickFromContextMenu(page, "Edit content...", el("markup viewer")));
      await session.step(39, "Then Edit dialog should be visible", () => shouldBe(page, el("Edit dialog"), "visible"));
      await session.step(40, "And input in Edit dialog should have the value \"editable probe\"", () => shouldHaveValue(page, el("input in Edit dialog"), "editable probe"));
      await session.step(41, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(42, "Then Edit dialog should be hidden", () => shouldBe(page, el("Edit dialog"), "hidden"));
      await session.step(43, "And the \"text\" reading of markup viewer should be \"editable probe\"", () => readingReads(page, "text", el("markup viewer"), "editable probe"));
      await session.step(44, "And \"content\" property of markup viewer should be \"editable probe\"", () => propertyShouldBe(page, "content", el("markup viewer"), "editable probe"));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The dialog opens on whatever the content is at the time, not on the default", async () => {
      await session.step(48, "When user sets \"content\" property of markup viewer to \"second probe\"", () => setProperty(page, "content", el("markup viewer"), "second probe"));
      await session.step(49, "And user picks \"Edit content...\" from the context menu of markup viewer", () => pickFromContextMenu(page, "Edit content...", el("markup viewer")));
      await session.step(50, "Then Edit dialog should be visible", () => shouldBe(page, el("Edit dialog"), "visible"));
      await session.step(51, "And input in Edit dialog should have the value \"second probe\"", () => shouldHaveValue(page, el("input in Edit dialog"), "second probe"));
      await session.step(52, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(53, "Then Edit dialog should be hidden", () => shouldBe(page, el("Edit dialog"), "hidden"));
      await session.step(54, "When user sets \"content\" property of markup viewer to \"editable probe\"", () => setProperty(page, "content", el("markup viewer"), "editable probe"));
      await session.step(55, "Then the \"text\" reading of markup viewer should be \"editable probe\"", () => readingReads(page, "text", el("markup viewer"), "editable probe"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title the property sets is the title the bar shows", async () => {
      await session.step(59, "Then title of markup viewer should not contain the text \"Patient card\"", () => shouldNotContainText(page, el("title of markup viewer"), "Patient card"));
      await session.step(60, "When user sets \"title\" property of markup viewer to \"Patient card\"", () => setProperty(page, "title", el("markup viewer"), "Patient card"));
      await session.step(61, "Then title of markup viewer should have text \"Patient card\"", () => shouldHaveText(page, el("title of markup viewer"), "Patient card"));
      await session.step(62, "When user sets \"title\" property of markup viewer to \"\"", () => setProperty(page, "title", el("markup viewer"), ""));
      await session.step(63, "Then title of markup viewer should not contain the text \"Patient card\"", () => shouldNotContainText(page, el("title of markup viewer"), "Patient card"));
      await session.step(64, "And the \"text\" reading of markup viewer should be \"editable probe\"", () => readingReads(page, "text", el("markup viewer"), "editable probe"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The title bar closes the markup viewer", async () => {
      await session.step(68, "When user clicks on close icon of markup viewer", () => clickOn(page, el("close icon of markup viewer")));
      await session.step(69, "Then markup viewer should be absent", () => shouldBe(page, el("markup viewer"), "absent"));
      await session.step(70, "And the open tableview should have 0 markup viewers", () => viewerCount(page, 0, "markup"));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The viewer's own menu opens anywhere on it, and the content fills it", async () => {
      await session.step(74, "Given user adds a markup viewer with:", () => addViewerWith(page, "markup", [["content","one line"]]));
      await session.step(76, "Then markup viewer should have a \"content\" area", () => hasArea(page, el("markup viewer"), "content"));
      await session.step(77, "And markup viewer should not have an \"empty space\" area", () => hasNoArea(page, el("markup viewer"), "empty space"));
      await session.step(78, "When user opens the viewer menu of markup viewer", () => openViewerMenu(page, el("markup viewer")));
      await session.step(79, "Then the open menu should list \"Edit content...\"", () => menuLists(page, "Edit content..."));
      await session.step(80, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(81, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
