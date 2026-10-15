---
feature: sequencetranslator
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: []
realizes: [sequencetranslator.app.oligo-toolkit]
realized_as: []
related_bugs: [GROK-16674, GROK-14385]
---

# SequenceTranslator — Oligo Pattern: edit, save, reload, delete

Checks the Oligo Pattern tool. A modification pattern is edited, and the example
translation follows the edits. The pattern is saved to the user's storage and
survives reopening the app. Deleting it removes the pattern chosen in the Load list.
The shipped example pattern cannot be saved over or deleted.

## Setup

- Open **Browse > Apps > Peptides > Oligo Toolkit > Oligo Pattern** (or open **Oligo
  Toolkit** and click the **PATTERN** tab).
- The pattern name used below is `BddPattern<N>`, where `<N>` is a number unique to
  the run. It must not exist yet.

## Scenarios

### Block A — Default example pattern

1. Look at the **Load** block.
   * Expected result: **Author** shows the current user's name followed by `(me)`.
     **Pattern** shows `<default example>`. The pattern picture and the **Translation
     example** block are shown, with examples for **Sense strand** and **Anti sense**.
     **Save** is disabled.
2. With `<default example>` selected, click the trash icon next to **Pattern**.
   * Expected result: warning balloon **Cannot delete example pattern**. No
     confirmation dialog.
3. In **Edit**, change **Sense strand length** to `22`, then click **Save**.
   * Expected result: **Save** becomes enabled after the change; clicking it shows the
     warning balloon **Cannot save default pattern**. The **Pattern** list is
     unchanged.

### Block B — Edit, save and reload

1. In **Edit**, set **Sense strand length** to `10` and switch **Anti sense strand**
   off.
   * Expected result: the picture shows a single strand. The **Anti sense** example
     disappears and the **Anti sense length** input is hidden. The **Sense strand**
     example input is 10 characters long.
2. Note the **Sense strand** example output. Click **Edit strands**. In the **Edit
   strands** dialog, switch **All PTO** on, change the modification of position 1
   from **RNA** to **2'-Fluoro**, and click **OK**.
   * Expected result: the **Sense strand** example output differs from the noted
     value.
3. Set **Pattern name** to `BddPattern<N>` and click **Save**.
   * Expected result: info balloon **Pattern BddPattern<N> saved**. The **Pattern**
     list in **Load** now contains `BddPattern<N>`.
4. Close the view and open the **Oligo Pattern** app again. In **Load**, choose
   `BddPattern<N>`.
   * Expected result: **Sense strand length** is `10`, **Anti sense strand** is off,
     and the **Sense strand** example output equals the one after step 2.

### Block C — Delete deletes the pattern chosen in Load (GROK-16674)

1. With `BddPattern<N>` loaded, change **Pattern name** in **Edit** to
   `OtherName<N>` (do not save).
2. Click the trash icon next to **Pattern** in **Load**.
   * Expected result: the **Delete pattern** dialog says **Are you sure you want to
     delete pattern BddPattern<N>?** — the name chosen in Load, not `OtherName<N>`.
3. Click **OK**.
   * Expected result: `BddPattern<N>` is no longer in the **Pattern** list. After
     reopening the app, it is still absent.

## Cleanup

- If a run stops before Block C step 3, open **Oligo Pattern**, choose
  `BddPattern<N>` in **Load**, click the trash icon and confirm with **OK**.

## Automation notes

- The pattern picture is an SVG drawing and is not checked; the **Translation
  example** text areas carry the observable effect of the edits.
- Block B step 2: the per-position modification inputs in **Edit strands** have no
  caption (only a position number next to them); the step needs a reading that finds
  the input for position 1 of the sense strand.

---
{
  "order": 3
}
