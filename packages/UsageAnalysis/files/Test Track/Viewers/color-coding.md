### 1. Setup
Open the **demog** table.

### 2. Apply color coding types
*(All steps on the same open table.)*

**2.1 AGE (numeric)**
- Apply Linear.
- Change to Conditional: `< 30` → green, `30–60` → yellow, `> 60` → red.

**2.2 SEX (string)**
- Apply Categorical with custom colors: `M` → `#3366CC`, `F` → `#CC6699`.

**2.3 CONTROL**
- Apply Categorical (default colors).

**2.4 STARTED (date)**
- Apply Linear with a custom 3-stop scheme: `#0000FF` → `#FFFFFF` → `#FF0000`.

Expected: each column shows its assigned type.

### 3. Disable and re-enable
**3.1** Turn off AGE, SEX, STARTED.
Expected: `getType()` = `"Off"` for each.

**3.2** Re-enable each column to its previous type.
Expected: custom colors from step 2 are preserved.

### 4. Pick Up / Apply coloring
**4.1** Create **Race_copy** column (copy of RACE values); apply Categorical.
**4.2** Apply RACE's color scheme to Race_copy.
**4.3** Apply STARTED's linear scheme to HEIGHT.

Expected: target columns have the same color type as their sources.

### 5. Linked color coding

**5.1 RACE linked to WEIGHT (background)**
- WEIGHT: apply Categorical.
- RACE: set type Linked, source = WEIGHT.

Expected: RACE background uses WEIGHT's categorical colors.

**5.2 HEIGHT linked to WEIGHT (text)**
- HEIGHT: set type Linked, Apply to: Text, source = WEIGHT.

Expected: HEIGHT text is colored using WEIGHT's colors.

**5.3 Source color changes propagate**
- Change WEIGHT to Linear → verify RACE and HEIGHT remain Linked.
- Change WEIGHT to Conditional → verify RACE and HEIGHT remain Linked.

Expected: linked columns automatically inherit WEIGHT's type changes.

**5.4 Linking chain (5 levels)**
- AGE (Linear) → SEX (Linked → AGE) → DIS_POP (Linked → SEX) → CONTROL (Linked → DIS_POP) → STARTED (Linked → CONTROL).

Expected: all 4 linked columns report type `"Linked"`.

### 6. Edit linear color scheme
**6.1** Apply Linear to AGE with a custom 3-stop scheme: `#1A237E` → `#F5F5F5` → `#B71C1C`.
**6.2** Invert the scheme.

Expected: gradient reverses (`#B71C1C` → `#F5F5F5` → `#1A237E`).

---
{
  "order": 1
}
