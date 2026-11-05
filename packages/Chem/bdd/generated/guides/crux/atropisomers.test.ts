/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/atropisomers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/datasets.js';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {cruxMenu, cruxOpenOnMolfile} from '../../../bindings/crux.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Atropisomers in Crux", () => {
  const session = feature(test, "features/guides/crux/atropisomers.feature", import.meta.url);
  test("Turn a biaryl's axis from P to M and back from its menu", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And the Crux sketcher is open on this molfile, showing R, S, E and Z labels:", () => cruxOpenOnMolfile(page, "\n     RDKit          2D\n\n 16 17  0  0  0  0  0  0  0  0999 V2000\n    0.6495    2.6250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6495    1.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9486    2.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.2476    1.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.2476    0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9486   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9486   -1.8750    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6495    0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6495   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9486    0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9486    1.8750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2476   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2476   -1.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9486   -2.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6495   -1.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6495   -2.6250    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  5  6  2  0\n  6  7  1  0\n  6  8  1  0\n  8  9  1  0\n  9 10  2  0\n 10 11  1  0\n 10 12  1  0\n 12 13  2  0\n 13 14  1  0\n 14 15  2  0\n 15 16  1  0\n  8  2  2  0\n  9 15  1  6\nM  END"));
    await session.step(57, "When user opens the Crux context menu on the \"bond 7\" area", () => cruxMenu(page, "bond 7"), undefined, "Right-click the axis, the bond between the two rings");
    await session.step(59, "Then Crux atropisomer P item should be checked", () => shouldBe(page, el("Crux atropisomer P item"), "checked"), undefined, "The axis is P: its menu has P checked");
    await session.step(61, "When user clicks on Crux atropisomer M item", () => clickOn(page, el("Crux atropisomer M item")), undefined, "Choose M");
    await session.step(63, "And user opens the Crux context menu on the \"bond 7\" area", () => cruxMenu(page, "bond 7"), undefined, "Right-click the axis again");
    await session.step(65, "Then Crux atropisomer M item should be checked", () => shouldBe(page, el("Crux atropisomer M item"), "checked"), undefined, "Now M is checked, and the label beside the axis reads (M)");
    await session.step(67, "When user clicks on Crux atropisomer P item", () => clickOn(page, el("Crux atropisomer P item")), undefined, "Choose P to turn it back");
    await session.step(69, "And user opens the Crux context menu on the \"bond 7\" area", () => cruxMenu(page, "bond 7"), undefined, "Right-click the axis once more");
    await session.step(71, "Then Crux atropisomer P item should be checked", () => shouldBe(page, el("Crux atropisomer P item"), "checked"), undefined, "P again");
  });
});
