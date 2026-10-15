/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/crux/dative-bonds.feature
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
import {cruxOpenOn} from '../../../bindings/crux.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, simpleModeOff, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dragAreaToArea, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Dative bonds in Crux", () => {
  const session = feature(test, "features/guides/crux/dative-bonds.feature", import.meta.url);
  test("Bind pyridine's nitrogen to a platinum atom with a dative bond", {tag: ["@guide", "@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And the Crux sketcher is open on \"c1ccncc1.[Pt]\"", () => cruxOpenOn(page, "c1ccncc1.[Pt]"));
    await session.step(15, "When user clicks on Crux dative bond tool", () => clickOn(page, el("Crux dative bond tool")), undefined, "Pick the dative bond tool");
    await session.step(17, "And user drags the \"atom 3\" area of Crux sketcher widget to the \"atom 6\" area", () => dragAreaToArea(page, "atom 3", el("Crux sketcher widget"), "atom 6"), undefined, "Drag from the nitrogen, the donor, to the platinum, the acceptor");
    await session.step(19, "And user clicks on Crux clean up button", () => clickOn(page, el("Crux clean up button")), undefined, "Clean Up tidies the drawing");
    await session.step(21, "Then the \"smiles\" reading of Crux sketcher widget should be \"[Pt]<-[n]1ccccc1\"", () => readingReads(page, "smiles", el("Crux sketcher widget"), "[Pt]<-[n]1ccccc1"), undefined, "An N→Pt dative bond: pyridine bound to platinum");
  });
});
