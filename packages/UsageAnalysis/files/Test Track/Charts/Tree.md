#### Tree Viewer

### 1. Collaborative Filtering

**Setup:**  
- Open **demog.csv**.
- Click **Add Viewer** → **Tree** to add the viewer.
- Set the **Hierarchy** to: **CONTROL**, **SEX**, and **RACE**.

---

**Steps:**

1. **Select branches** in the tree (use **Shift + Click**):  
   - **All → false → F → Asian**  
   - **All → false → F → Black**  
   - **All → false → M → Asian**

2. In the **Filter panel**, select:  
   - **CONTROL = true**

**Expected**: Filtered count = 0

---

3. **Modify selection** in the tree:  
   - Add **All → true → F → Black** (use **Shift + Click** to multi-select).

**Expected**: Filtered count = 2

---

4. In the **Filter Panel**, **clear** the filter for **CONTROL = true**.

**Expected**: Filtered count = 176

---
{
  "order": 30,
  "datasets": ["System:DemoFiles/demog.csv"]
}