/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/15-annotate-charts/07-add-a-note.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, hoverOver, shouldBe, shouldContainText, typeInto, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, noErrors, pickFromContextMenu, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Annotations", () => {
  const session = feature(test, "features/howto/15-annotate-charts/07-add-a-note.feature", import.meta.url);
  test("Add a note that explains a chart", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(18, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(19, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(20, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(23, "When the video says \"A Markup viewer is a note next to your charts\"", () => videoSays(page, "A Markup viewer is a note next to your charts"));
    await session.step(25, "And user clicks on \"Add viewer\" icon in toolbar", () => clickOn(page, el("\"Add viewer\" icon in toolbar")), undefined, "Add viewer, then search for Markup");
    await session.step(26, "And user types \"markup\" into viewer gallery search in \"Add Viewer\" dialog", () => typeInto(page, "markup", el("viewer gallery search in \"Add Viewer\" dialog")));
    await session.step(27, "And user clicks on first \"Markup\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Markup\" card in \"Add Viewer\" dialog")));
    await session.step(28, "Then markup viewer should be visible", () => shouldBe(page, el("markup viewer"), "visible"));
    await session.step(30, "When user picks \"Edit content...\" from the context menu of markup viewer", () => pickFromContextMenu(page, "Edit content...", el("markup viewer")), undefined, "Right-click it: Edit content...");
    await session.step(31, "Then Edit dialog should be visible", () => shouldBe(page, el("Edit dialog"), "visible"));
    await session.step(32, "When user enters \"## Patient ${USUBJID}: age ${AGE}, ${WEIGHT} kg, ${HEIGHT} cm\" into input in Edit dialog", () => enterInto(page, "## Patient ${USUBJID}: age ${AGE}, ${WEIGHT} kg, ${HEIGHT} cm", el("input in Edit dialog")));
    await session.step(34, "And user hovers over input in Edit dialog", () => hoverOver(page, el("input in Edit dialog")), undefined, "Markdown text; ${COLUMN} shows the current row's value");
    await session.step(35, "And user clicks on \"OK\" button in Edit dialog", () => clickOn(page, el("\"OK\" button in Edit dialog")));
    await session.step(36, "Then Edit dialog should be absent", () => shouldBe(page, el("Edit dialog"), "absent"));
    await session.step(37, "And \"content\" property of markup viewer should be \"## Patient ${USUBJID}: age ${AGE}, ${WEIGHT} kg, ${HEIGHT} cm\"", () => propertyShouldBe(page, "content", el("markup viewer"), "## Patient ${USUBJID}: age ${AGE}, ${WEIGHT} kg, ${HEIGHT} cm"));
    await session.step(39, "When user clicks on the \"marker of row 11\" area of scatter plot viewer", () => clickArea(page, "marker of row 11", el("scatter plot viewer")), undefined, "Click a point: the note shows that patient");
    await session.step(41, "Then row 11 should be current", () => currentRowIs(page, 11), undefined, "Age 46, 118.90 kg, 168.974 cm: the point you clicked");
    await session.step(42, "And markup viewer should contain text \"Patient X0273T21000900001: age 46, 118.90 kg, 168.974 cm\"", () => shouldContainText(page, el("markup viewer"), "Patient X0273T21000900001: age 46, 118.90 kg, 168.974 cm"));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
});
