# Browse — Apps > Demo (manual test cases)

Мануальные тест-кейсы для просмотра **всех демо-приложений** в `Browse > Apps > Demo`.
Каждое демо открывается **кликом по узлу дерева Browse (по названию)**, после чего во вью
выполняется минимальное взаимодействие. Главная проверка — **отсутствие любых ошибок**.

Автоматизация: `demo_apps.test.ts` (этот же каталог). Кейсы структурированы под прямой перенос
в Playwright — матрица демо ниже 1:1 совпадает с массивом `DEMOS` в коде.

---

## 1. Окружение и общие предусловия

- **Инстанс:** `https://dev.datagrok.ai/` под залогиненным пользователем с правами на демо-данные.
- **Источник дерева:** функции с `meta.demoPath` (формат `Категория | Подкатегория | Демо`),
  снято живьём с dev — **62 демо** в 6 категориях. Список — раздел 4.
- **Открыто:** боковая панель **Browse** (иконка на Sidebar), узел **Apps** развёрнут.
- Навигация строго **по названию** узла, не по индексу/порядку.

## 2. Как открывается демо (поведение платформы)

- Клик по листу демо вызывает функцию демо (`func.apply()`), которая открывает новый view и
  делает его текущим. **Имя view = последний сегмент пути** (напр. `Scatter Plot`).
- Рядом обычно открывается вторичный data-view (`Table` / `Grid`) с исходными данными.
- В шапке view — breadcrumbs `Home / Demo / … / <Демо>`.
- Открытие следующего демо закрывает предыдущие demo-view (они помечены `temp.demoApp`).
- Часть демо — **дашборды** (`isDemoDashboard`): открываются как сохранённый layout.
- **Тяжёлые** демо (Compute, Docking, Admetica, Retrosynthesis, Databases, Map…) запускают
  серверные/докер-вычисления — открываются дольше, могут зависеть от доступности бэкенда на dev.

## 3. Общий тест-кейс (применяется к каждому демо)

**ID:** `Browse-DemoApps-<Категория>-<Демо>`
**Тип:** smoke (для всех); regression — для демо с `ref` (раздел 4).

**Предусловия:** раздел 1.

**Шаги (действие → ожидаемая реакция):**

1. В дереве Browse раскрыть `Apps` → реакция: появляются дочерние узлы, среди них `Demo`.
2. Раскрыть `Apps > Demo` → появляются категории (Cheminformatics, Bioinformatics, Data Access,
   Visualization, Compute, Curves).
3. Раскрыть категорию (и подкатегорию, если есть, напр. `Visualization > General`) **по названию**
   → появляются листья-демо.
4. Кликнуть лист демо **по названию** → во вью открывается view с заголовком = имя демо;
   в шапке — breadcrumbs `Home / Demo / … / <Демо>`.
5. **Поклацать во вью:** кликнуть по телу view; если есть таблица — кликнуть по ячейке grid
   (меняется текущая строка); навести курсор на основной viewer → реакция: выделение/тултип,
   перерисовка без сбоев.
6. Открыть **Context Panel** (F4) и нажать **Expand all** → все панели свойств отрисовываются.

**Итоговая проверка (для всех демо):**

- ✅ view с именем демо открылся;
- ✅ нет uncaught JS-ошибок (`pageerror`);
- ✅ нет `console.error` (кроме явно занесённого в игнор-лист шума — см. `helpers.ts`);
- ✅ нет error-баллунов;
- ✅ при раскрытии Context Panel ни одна панель не падает (ошибка панели всплыла бы в
  console/pageerror).

**Постусловия / очистка:** вернуться `Home` (закрывает demo-preview); ничего на сервере не
создаётся — серверная очистка не нужна.

## 4. Матрица демо (62)

Колонки: **Path** — путь под `Apps > Demo`; **View** — имя открытого view; **Tags** —
`dashboard` (layout), `heavy` (сервер/докер), `skip` (платформа метит `demoSkip`); **Ref** —
тикет/флаг.

### Cheminformatics (10)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Cheminformatics / Med Chem | Med Chem | | |
| Cheminformatics / Chemical Space | Chemical Space | skip | GROK-14320 |
| Cheminformatics / Molecule Activity Cliffs | Molecule Activity Cliffs | skip, dashboard | GROK-14320 |
| Cheminformatics / R-Group Analysis | R-Group Analysis | skip, dashboard | GROK-14320 |
| Cheminformatics / Matched Molecular Pairs | Matched Molecular Pairs | | |
| Cheminformatics / Similarity & Diversity Search | Similarity & Diversity Search | | |
| Cheminformatics / Scaffold Tree | Scaffold Tree | | |
| Cheminformatics / Database Queries | Database Queries | heavy | |
| Cheminformatics / Admetica | Admetica | heavy | |
| Cheminformatics / Retrosynthesis | Retrosynthesis | heavy | |

### Bioinformatics (10)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Bioinformatics / Peptide SAR | Peptide SAR | dashboard, heavy | |
| Bioinformatics / Sequence Activity Cliffs | Sequence Activity Cliffs | heavy | |
| Bioinformatics / siRNA | siRNA | skip, heavy | demoSkip |
| Bioinformatics / Antibodies | Antibodies | heavy | |
| Bioinformatics / Sequence Space | Sequence Space | dashboard, heavy | |
| Bioinformatics / Similarity, Diversity | Similarity, Diversity | heavy | |
| Bioinformatics / Atomic Level | Atomic Level | skip, heavy | demoSkip |
| Bioinformatics / Docking | Docking | heavy | |
| Bioinformatics / Docking Conformations | Docking Conformations | heavy | |
| Bioinformatics / Proteins | Proteins | heavy | |

### Data Access (3)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Data Access / Table Linking | Table Linking | | |
| Data Access / Files | Files | | |
| Data Access / Databases | Databases | heavy | |

### Visualization (33)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Visualization / Data Flow and Hierarchy / Network Diagram | Network Diagram | | |
| Visualization / Data Flow and Hierarchy / Tree | Tree | | |
| Visualization / Data Flow and Hierarchy / Tree Map | Tree Map | heavy | |
| Visualization / Data Separation / Trellis Plot | Trellis Plot | | |
| Visualization / Data Separation / Matrix Plot | Matrix Plot | | |
| Visualization / General / Scatter Plot | Scatter Plot | | |
| Visualization / General / Bar Chart | Bar Chart | | |
| Visualization / General / Line Chart | Line Chart | | |
| Visualization / General / Histogram | Histogram | | |
| Visualization / General / Pie Chart | Pie Chart | | |
| Visualization / General / 3D Scatter Plot | 3D Scatter Plot | | |
| Visualization / General / Tile Viewer | Tile Viewer | | |
| Visualization / General / Density Plot | Density Plot | | |
| Visualization / General / Filters | Filters | | |
| Visualization / General / Heatmap | Heatmap | | |
| Visualization / General / Markup | Markup | | |
| Visualization / General / Radar | Radar | | |
| Visualization / General / Sunburst | Sunburst | | |
| Visualization / General / Chord | Chord | | |
| Visualization / General / Sankey | Sankey | | |
| Visualization / General / Surface Plot | Surface Plot | | |
| Visualization / General / Timelines | Timelines | | |
| Visualization / General / Word Cloud | Word Cloud | | |
| Visualization / General / Data Annotations | Data Annotations | | |
| Visualization / Geographical / Map | Map | heavy | |
| Visualization / Input and Edit / Grid | Grid | | |
| Visualization / Input and Edit / Form | Form | | |
| Visualization / Statistical / Box Plot | Box Plot | | |
| Visualization / Statistical / Correlation Plot | Correlation Plot | | |
| Visualization / Statistical / PC Plot | PC Plot | | |
| Visualization / Statistical / Pivot Table | Pivot Table | | |
| Visualization / Statistical / Statistics | Statistics | | |
| Visualization / Time and Date / Calendar | Calendar | | |

### Compute (4)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Compute / Diff Studio | Diff Studio | heavy | |
| Compute / Multivariate Analysis | Multivariate Analysis | heavy | |
| Compute / PK-PD Modeling | PK-PD Modeling | heavy | |
| Compute / Bioreactor | Bioreactor | heavy | |

### Curves (2)

| Path | View | Tags | Ref |
|------|------|------|-----|
| Curves / Curve Fitting | Curve Fitting | | |
| Curves / Assay Curves | Assay Curves | | |

## 5. Результаты первого прогона на dev (триаж)

Прогон всех 62 демо на `dev.datagrok.ai` (2026-06-04): **59 passed, 1 known-fail, 2 поправлено
в тесте**.

- ✅ **59 демо** открываются и кликаются без ошибок — включая все «skip» (GROK-14320:
  Chemical Space, Molecule Activity Cliffs, R-Group Analysis; `demoSkip`: siRNA, Atomic Level)
  и все тяжёлые Compute/Docking/Admetica/Retrosynthesis. Т.е. `demoSkip` относится к
  demo-runner'у, а не к Browse.
- 🐞 **Bioinformatics / Similarity, Diversity** — «случайные клики» по вью (шаг `pokeView`)
  бросают `NullError: method not found: 'gdV' on null`. Это **регрессия
  [GROK-18050](https://reddata.atlassian.net/browse/GROK-18050)** (тикет в статусе Done).
  Тест **обычный, без подавления** — падает по-настоящему, когда баг срабатывает (он
  перемежающийся, т.к. клик по фиксированным координатам не всегда попадает на битый элемент).
  Тикет привязан через `test.info().annotations`. Подробности — `KNOWN_BUGS.md`.
- 🔧 **Med Chem** — `console.error` был от CORS-блокировки внешней картинки Wikipedia в
  описании демо (не дефект платформы) → добавлен в игнор-лист `helpers.ts`.
- 🔧 **Database Queries** — открывает форму SQL-запроса, а не data-view; помечено `opensForm`.
  При авто-прогоне на dev демо-БД Chembl (`db.datagrok.ai:54325`) была недоступна
  (`Connection refused`) — это **инфраструктура dev**, не дефект Browse/демо. Такая ошибка
  бэкенда толерируется через `envIgnore` (любая другая ошибка по-прежнему валит тест). Сама
  навигация и открытие формы работают.

## 6. Заметки и известные проблемы

- **`skip` (`demoSkip`)** — платформа исключает эти демо из автоматического demo-runner'а как
  нестабильные. В наборе они **не пропускаются**: цель — поймать реальные ошибки. Те, что
  стабильно падают на dev, фиксируются в `KNOWN_BUGS.md` и помечаются `test.fail()` (как
  ModelHub-кейсы), чтобы фикс на сервере сразу всплыл как «unexpected pass».
- **GROK-14320** — Chemical Space, Molecule Activity Cliffs, R-Group Analysis (Cheminformatics).
- **`heavy`** — может падать не из-за UI, а из-за недоступности докер-бэкенда на dev
  (Admetica, Retrosynthesis, Docking, Boltz и т.п.). Такие падения трактуем отдельно от
  ошибок Browse/рендера.
- **Plates** — в источниках есть `Plates | Assay Plates`, но на dev в дереве Demo **отсутствует**
  (не задеплоено). В матрицу не включено.
- Список снят автоматически; при изменении набора демо на сервере перегенерировать матрицу и
  массив `DEMOS`.
