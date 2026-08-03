# Datagrok Browse — контекст проекта (handoff)

Этот документ — сопровождение к `browse_manual_tests.md`. Содержит весь контекст, который нужен, чтобы продолжить работу в новом чате: цель, окружение, решения по скоупу, источники, выжимку по устройству Browse и план следующих шагов.

---

## 1. Цель и текущий этап

- Пишем **мануальные тест-кейсы для Datagrok Browse**, покрытие через UI.
- Конечная цель — **автотесты на Playwright**. Мануалка специально структурирована под лёгкий перенос в код.
- Сейчас этап: готова мануалка (`browse_manual_tests.md`, 89 кейсов, 17 секций). Дальше — снять селекторы на dev и начать перенос в Playwright.

## 2. Окружение

- **Целевой инстанс:** https://dev.datagrok.ai/
- Тесты под залогиненного пользователя с правами на демо-данные.
- На dev есть **демо-файлы** (Files > Demo), **демо-БД** (Databases), **демо-дашборды** (Dashboards) — предусловия опираются на них.

## 3. Формат тест-кейсов

- Markdown, текстовое описание (не таблицы).
- Каждый кейс: **Цель → Тип → Предусловия → Шаги (действие + ожидаемая реакция) → Итоговая проверка → Постусловия/очистка**.
- **ID:** `Browse-<Область>-<NN>` — формат совместим с Datagrok Test Track (ср. `Scripts-Edit-6`).
- **Тип:** `smoke` / `functional` / `regression` / `negative` — пойдёт в теги Playwright (`@smoke` и т.д.).
- Шаги атомарные: действие → `await`-вызов, ожидание → `expect`. Предусловия → `beforeEach`/фикстура, постусловия → `afterEach`/teardown.

## 4. Решения по скоупу (согласовано)

- **Покрываем всё дерево Browse**: My stuff, Spaces, Apps, Files, Dashboards, Databases, Platform + сквозные (навигация, дерево, view-режимы, фильтры/поиск, routing).
- **Spaces** — в Browse только **смоук** (клик по ноде без ошибок в консоли). Полное покрытие Spaces (rename, перемещение, наследование прав, потеря данных) — в **отдельном наборе тестов по Spaces**, не здесь.
- **Model Hub / Compute** — только **смоук** в Browse (точки исторических падений).
- **Global search, Context Panel controls, Split view / viewing modes** — добавлены как базовые секции (15–17).

## 5. Где живут тест-кейсы в Datagrok

- В Datagrok есть встроенный инструмент **Test Track** (DevTools → Test Track) для мануальных тест-кейсов. ref: GROK-14402 (DevTools: Test Track: Manual testing).
- Команда движется к замене тикетов полноценными мануальными кейсами. ref: GROK-15857 (Replace tickets in the test track with manual test cases) — наша работа ложится в этот подход.
- Кейсы в Test Track именуются по шаблону `<Область>-<Действие>-<номер>` — мы используем совместимый формат.

## 6. Источники

- Документация: `https://datagrok.ai/help/datagrok/navigation/` (раздел Navigation / Browse).
- Jira: инстанс `reddata.atlassian.net`, проект **GROK** (Datagrok). cloudId: `42e2928b-7dc3-47f5-8ab7-52961685645a`.
  - Доступ к Jira: работает JQL-поиск (`searchJiraIssuesUsingJql`) и чтение issue; полнотекстовый Rovo `search` отдавал 403 («app is not installed») — использовать JQL.

## 7. Баги-источники (привязка regression-кейсов)

Список багов Browse, на которые ссылаются кейсы (поле `ref:`). Полезны и как регрессия, и как карта «где исторически ломалось».

- **GROK-16261** — browse tree должно запоминать состояние expand/collapse (Tree-05).
- **GROK-19802** — вложенные узлы «залипают» раскрытыми при сворачивании бокового меню (Tree-06).
- **GROK-17922** — клик по объекту без прав: должна быть заглушка, не падение (Tree-09, negative).
- **GROK-19848** — Add to Favorites из MyFiles (MyStuff-05, Fav).
- **GROK-20032** — список Apps грузился ~40 сек (Apps-02).
- **GROK-19638** — ошибки в тултипе/details приложения (Apps-04).
- **GROK-19562** — открытие .h5 файла (Files-04).
- **GROK-19844** — новый файл в S3/Spaces появлялся только после ручного refresh (Files-05).
- **GROK-19847** — иконка share-папки (Files-06).
- **GROK-19662** — новый дашборд появлялся только после refresh (Dash-03).
- **GROK-19934** — дублирующийся контент в Context Panel при смене дашборда (Dash-04).
- **GROK-16857** — Postgres > CHEMBL > Browse > Summary падал (DB-06).
- **GROK-19691** — пустые/нерелевантные свойства в панели фильтров (Filter-01).
- **GROK-19689** — фильтр CreatedRecently для Files (Filter-02).
- **GROK-19688** — фильтр UsedByMe для Apps/Dockers (Filter-03).
- **GROK-19692** — фильтр Users по группе, ошибки/«балуны» (Filter-04).
- **GROK-19690** — совместная работа фильтра и поиска (Filter-05).
- **GROK-19703** — фильтр Models по Namespace (Filter-06).
- **Model Hub:** GROK-19965 (двойной клик падал, right-click→Run работал), GROK-19740 (клик/hover по модели), GROK-19628 (раскрытие Uncategorized), GROK-17896 / GROK-17664 (открытие Model Catalog из Apps), GROK-17770 (открытие после сохранения).

> Список не исчерпывающий. При расширении можно искать в Jira: `project = GROK AND summary ~ "Browse" ORDER BY created DESC`. Также смотреть активные баги по Spaces (GROK-19843/54/57/58, 19845) — но они идут в отдельный набор Spaces-тестов.

## 8. Выжимка по устройству Browse (из доки)

**Что такое Browse:** единая точка входа ко всему — projects, queries, files, connections, settings, apps. Открывается иконкой на Sidebar.

**Дерево Browse (секции верхнего уровня):**
- **My stuff** — Recent, Favorites, Shared with me, My Dashboards, My Files и др. Объект без явного проекта попадает сюда.
- **Spaces** — общие пространства (по правам).
- **Apps** — установленные приложения.
- **Files** — My files / Spaces / App Data / Demo.
- **Dashboards** — свои + расшаренные.
- **Databases** — connections → tables → columns → queries; есть Schema Browser.
- **Platform** — конфигурация инстанса (плагины, функции, access control); по правам.

**Top Menu панели Browse:** Home, Import file, Import text, Refresh tree, Collapse all, Locate current object.

**Browsing mode vs Persistent view:**
- Одиночный клик → **unpinned** view (заменяется следующим кликом).
- Двойной клик / любое редактирование / pin → **persistent** view (не заменяется).
- На Sidebar однотипные view группируются с бейджем-счётчиком.

**Управление деревом:**
- Клик / двойной клик / правый клик (контекстное меню).
- Клавиатура: ↑/↓ — навигация, →/← — expand/collapse, Enter/Space — открыть элемент.
- Drag-and-drop — реорганизация дерева и перетаскивание сущностей в панели.

**Context Panel** (правый экран): тоггл **F4**; в шапке — Back/Forward, Clone and detach, Collapse all / Expand all, Favorites (звезда). Контент следует за выбранным объектом.

**Split & viewing modes:** Hamburger → Split right / Split down; режимы со Status Bar — Tabbed / Simple / Presentation.

**Global search:** строка на Home Page; поиск по метаданным, тегам (`#tag`), идентификаторам, есть NLP.

**Favorites:** доступно только для **entities** (не для отдельных файлов и не для значений ячеек). Добавление: контекстное меню, звезда в Context Panel, drag-and-drop в Favorites Panel.

**Routing:** URL отражает состояние — дашборды, файлы (`/files/{namespace}.{share}/{path}`), параметризованные queries (`/q/{ns}.{conn}.{query}?param=value`).

## 9. Заметки для автоматизации (сводка)

- **baseURL** = `https://dev.datagrok.ai/`; авторизация через `storageState`/глобальный setup.
- **assertNoErrors(page)** — общий хелпер: `page.on('pageerror')` + `page.on('console')` (фильтр error) + проверка отсутствия datagrok error-балуна. Звать в финале всех «без ошибок»-кейсов.
- **Создание данных для предусловий** (Nav-06, Files-05, MyStuff-04, ModelHub-05) — через API/второй контекст, не руками в UI.
- **Очистка** — имена с префиксом `qa_autotest_<timestamp>` + общий teardown, удаляющий всё по префиксу (страховка от упавших тестов).
- **Асинхронность** — ждать появления узлов/элементов, не `sleep`.
- **Тайминги** — порог Apps-02 (5 c) вынести в константу конфига.
- **Селекторы пока не сняты.** На первом прогоне на dev снять стабильные `data-*`/роли для: иконок Browse/Favorites на Sidebar; узлов дерева и стрелок expand/collapse; Top Menu; контекстного меню; звезды Favorites и контролов Context Panel; табов и кнопки pin; Hamburger (Split); переключателя режимов на Status Bar; строки глобального поиска; панели фильтров и поля поиска.

## 10. План дальнейших шагов

1. Пройти **P1 smoke**-кейсы руками на dev → снять реальные селекторы → добавить блок «Селекторы» в каждый кейс.
2. Согласовать набор с командой; при необходимости завести в Test Track.
3. Перенос в Playwright по очереди: сперва `@smoke` (P1), затем `@regression` (баговые кейсы), затем остальное.

## 11. Открытые вопросы на потом

- Точные имена демо-сущностей на dev (файл/БД/дашборд), которые стабильно есть, — зафиксировать в фикстурах.
- Порог времени для Apps-02 (сейчас 5 c как предположение) — подтвердить с командой.
- Точный состав контекстного меню для разных типов сущностей (Tree-07) — снять на первом прогоне.
- Нужна ли отдельная секция Spaces целиком (сейчас договорились — отдельный набор тестов вне этого файла).
