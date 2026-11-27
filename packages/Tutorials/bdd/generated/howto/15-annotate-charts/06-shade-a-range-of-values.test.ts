/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/15-annotate-charts/06-shade-a-range-of-values.feature
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
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hasArea, hoverArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Annotations", () => {
  const session = feature(test, "features/howto/15-annotate-charts/06-shade-a-range-of-values.feature", import.meta.url);
  test("Shade a range of values on a chart", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
    await session.step(24, "And histogram viewer takes 65% of the space it shares with the grid", () => viewerShare(page, el("histogram viewer"), 65));
    await session.step(25, "When the video says \"Shade a range of values to point it out\"", () => videoSays(page, "Shade a range of values to point it out"));
    await session.step(27, "And user picks \"Annotations > Add Region\" from the context menu of the \"x axis\" area of histogram viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Region", "x axis", el("histogram viewer")), undefined, "Right-click the X axis: Annotations > Add Region");
    await session.step(28, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(30, "When user enters \"Middle half\" into \"Title\" text input in \"Formula Lines\" dialog", () => enterInto(page, "Middle half", el("\"Title\" text input in \"Formula Lines\" dialog")), undefined, "It spans the first to the third quartile: the middle half");
    await session.step(31, "And user clicks on OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(32, "Then histogram viewer should have a 'region \"Middle half\"' area", () => hasArea(page, el("histogram viewer"), "region \"Middle half\""));
    await session.step(34, "When user picks \"Annotations > Add Region\" from the context menu of the \"x axis\" area of histogram viewer", () => pickFromAreaContextMenu(page, "Annotations > Add Region", "x axis", el("histogram viewer")), undefined, "Add another, then type its bounds: AGE from 65 to 89");
    await session.step(35, "Then \"Formula Lines\" dialog should be visible", () => shouldBe(page, el("\"Formula Lines\" dialog"), "visible"));
    await session.step(36, "When user enters \"${AGE} = 65\" into \"Formula 1\" text area in \"Formula Lines\" dialog", () => enterInto(page, "${AGE} = 65", el("\"Formula 1\" text area in \"Formula Lines\" dialog")));
    await session.step(37, "And user enters \"${AGE} = 89\" into \"Formula 2\" text area in \"Formula Lines\" dialog", () => enterInto(page, "${AGE} = 89", el("\"Formula 2\" text area in \"Formula Lines\" dialog")));
    await session.step(38, "And user enters \"65 and over\" into \"Title\" text input in \"Formula Lines\" dialog", () => enterInto(page, "65 and over", el("\"Title\" text input in \"Formula Lines\" dialog")));
    await session.step(39, "And user clicks on OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(40, "Then histogram viewer should have a 'region \"65 and over\"' area", () => hasArea(page, el("histogram viewer"), "region \"65 and over\""));
    await session.step(41, "And histogram viewer should have a 'region \"Middle half\"' area", () => hasArea(page, el("histogram viewer"), "region \"Middle half\""));
    await session.step(43, "When user hovers over the 'region \"65 and over\"' area of histogram viewer", () => hoverArea(page, "region \"65 and over\"", el("histogram viewer")), undefined, "Point at a region: its name and how many rows fall in it");
    await session.step(44, "Then tooltip should contain text \"85 rows\"", () => shouldContainText(page, el("tooltip"), "85 rows"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
  });
});
