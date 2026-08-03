# Home page Widgets — мануальные тест-кейсы

Покрытие **виджетов на Home page** (лендинг-вью Datagrok): рендер, управление (close/Customize),
persistence, поиск, навигация и функционал отдельных виджетов (Spotlight, Community, Usage/Reports).

Контекст, окружение, селекторы и подводные камни — в `widgets_tests_context.md`.
Кейсы укрупнены и структурированы под прямой перенос в Playwright (один кейс = один сквозной поток).

---

## 1. Окружение и общие предусловия

- **Инстанс:** `https://dev.datagrok.ai/` под залогиненным пользователем (`opavlenko+playwright`, имеет роль admin).
- Открыт **Home Page** (логин или иконка **Home**). Видны виджеты **Spotlight**, **Community**, **Usage**, **Reports**.
- Главная сквозная проверка для всех кейсов — **отсутствие ошибок** (`pageerror`, `console.error`, error-балуны).
- **Очистка:** кейсы, меняющие видимость виджетов, восстанавливают исходное состояние (все виджеты `ignored:false`).

---

## 2. Display — рендер Home page и виджетов

### Widgets-Display-01 — Home page рендерится со всеми виджетами
**Тип:** smoke
**Предусловия:** раздел 1.
**Шаги:**
1. Залогиниться / открыть Home page → реакция: открывается вью `.power-pack-welcome-view`, появляется `.power-pack-widgets-host` с дочерними `.power-pack-widget-host`.
2. Осмотреть отрисованные виджеты (presence, порядок, заголовки, контент).
**Итоговая проверка:**
- ✅ присутствует строка поиска (placeholder начинается с «Search everywhere»);
- ✅ видны виджеты **Spotlight**, **Community**, **Usage**, **Reports** (`[widget-title="..."]`);
- ✅ **порядок:** Spotlight (order -1) идёт раньше Community (order 6) в DOM-порядке;
- ✅ **заголовки:** Community показывает `.d4-dialog-title` = «Community»; Spotlight — без заголовка (`showName=false`, `.d4-dialog-header-hidden`);
- ✅ контент каждого виджета (`.power-pack-widget-content`) непустой;
- ✅ нет ошибок (упавший виджет логирует `console.error` и удаляет свой хост — этого быть не должно).
**Постусловия:** —

---

## 3. Controls — управление виджетом (hover + close)

### Widgets-Controls-01 — Hover показывает close-иконку, клик убирает виджет
**Тип:** functional
**Предусловия:** Home page открыт, Community видим.
**Шаги:**
1. Навести курсор на виджет **Community** → реакция: становится видна close-иконка `i.grok-font-icon-close` (тултип «Remove»).
2. Кликнуть close-иконку → реакция: виджет исчезает из `.power-pack-widgets-host`.
**Итоговая проверка:**
- ✅ при hover close-иконка видима, без hover — скрыта;
- ✅ после клика `[widget-title="Community"]` отсутствует в DOM, остальные виджеты остаются.
**Постусловия / очистка:** вернуть Community через «Customize widgets...» → `ignored:false`.
> Авто-заметка: синтетический `mouseenter` не включает видимость иконки — использовать Playwright `.hover()`; сам клик по close работает и без hover (фолбэк `click({force:true})`).

### Widgets-Controls-02 — Закрытый виджет остаётся скрытым после reload (persistence)
**Тип:** regression
**Предусловия:** Home page открыт, Community видим.
**Шаги:**
1. Закрыть **Community** (close-иконка).
2. **Дождаться** сохранения настроек (async — см. авто-заметку).
3. Перезагрузить страницу (F5).
**Итоговая проверка:**
- ✅ после reload **Community отсутствует**;
- ✅ `grok.userSettings.get('widgets')` → `communityWidget` = `{"ignored":true}`.
**Постусловия / очистка:** вернуть Community через Customize → `ignored:false`, перепроверить reload.
> Авто-заметка: без паузы перед reload кейс флакает (сохранение `ignored:true` не долетает до сервера). Снято вживую 2026-06-08.

---

## 4. Customize — «Customize widgets...» (Context Panel)

### Widgets-Customize-01 — Форма настроек: скрытие и показ виджета
**Тип:** functional
**Предусловия:** Home page открыт, Community видим.
**Шаги:**
1. Кликнуть ссылку **«Customize widgets...»** под виджетами → реакция: в Context Panel открывается форма bool-тогглов (по тогглу на виджет, label = friendlyName; видимые отмечены).
2. Снять галку **Community** → реакция: виджет исчезает с Home page.
3. Поставить галку **Community** обратно → реакция: виджет снова появляется.
**Итоговая проверка:**
- ✅ в форме есть тогглы минимум **Spotlight, Community, Usage, Reports**;
- ✅ после снятия галки `[widget-title="Community"]` отсутствует, `communityWidget` = `{"ignored":true}`;
- ✅ после возврата галки виджет виден, `communityWidget` = `{"ignored":false}`.
**Постусловия:** состояние восстановлено.
> Авто-заметка: тот же тоггл возвращает виджет, ранее закрытый через close-иконку (`Widgets-Controls-01`) — отдельный кейс не нужен, механизм идентичен.

### Widgets-Customize-02 — Настройка видимости сохраняется после reload
**Тип:** regression
**Предусловия:** Home page открыт.
**Шаги:**
1. В Customize снять галку **Community**, дождаться сохранения.
2. Перезагрузить страницу.
**Итоговая проверка:**
- ✅ после reload Community скрыт;
- ✅ форма Customize по-прежнему показывает Community unchecked.
**Постусловия / очистка:** вернуть галку → `ignored:false`.

---

## 5. Search — взаимодействие поиска и виджетов

### Widgets-Search-01 — Ввод запроса прячет виджеты, очистка возвращает
**Тип:** functional
**Предусловия:** Home page открыт, панель виджетов видна.
**Шаги:**
1. В строку поиска ввести `aspirin` → (debounce ~500 мс) реакция: `.power-pack-widgets-panel` скрывается, `.power-pack-search-host` показывается с результатами.
2. Очистить строку поиска → реакция: панель виджетов снова видна, результаты скрыты.
**Итоговая проверка:**
- ✅ при вводе: панель виджетов скрыта, хост результатов виден, path = `search?q=aspirin`;
- ✅ при очистке: панель виджетов видна, path = `search` (без `?q=`).
**Постусловия:** —

---

## 6. Navigation — навигация на Home page

### Widgets-Nav-01 — Иконка Home возвращает на Home, набор виджетов стабилен
**Тип:** smoke
**Предусловия:** Открыт другой вью (например, таблица/Browse-нода).
**Шаги:**
1. Кликнуть иконку **Home** (`[name="icon-home"]`) → реакция: открывается Home page с виджетами.
2. Перезагрузить страницу.
**Итоговая проверка:**
- ✅ присутствует `.power-pack-welcome-view` с виджетами;
- ✅ набор видимых виджетов идентичен до/после reload (Spotlight, Community, Usage, Reports).
**Постусловия:** —

---

## 7. Spotlight — мульти-таб виджет

### Widgets-Spotlight-01 — Рендер с табами и «Demo of the day»
**Тип:** smoke
**Предусловия:** Home page открыт, Spotlight виден.
**Шаги:**
1. Осмотреть **Spotlight** → видны табы **Workspace, Spotlight, Favorites, Notifications, My Activity, Learn** и ссылка **«💡 Demo of the day: <name>»** внизу.
**Итоговая проверка:**
- ✅ присутствуют все 6 табов (`.d4-tab-header`);
- ✅ таб **Notifications** показывает бейдж-счётчик непрочитанных (если есть);
- ✅ ссылка «Demo of the day» присутствует и ведёт на демо при клике.
**Постусловия:** —

### Widgets-Spotlight-02 — Переключение табов (вкл. под-вкладки Learn)
**Тип:** functional
**Предусловия:** Spotlight виден.
**Шаги:**
1. Кликнуть **Workspace** → закреплённые (pinned) сущности, подсказка «Select a pinned item...».
2. Кликнуть **Spotlight** → Shared with me / Recent сущности.
3. Кликнуть **Favorites** → избранное.
4. Кликнуть **My Activity** → активность пользователя.
5. Кликнуть **Learn** → под-вкладки **VIDEO / WIKI / DEMO / TUTORIALS** и темы (Meetings, Develop, Cheminformatics, Visualize, Explore); переключить под-вкладки.
**Итоговая проверка:**
- ✅ каждый таб переключает контент `.power-pack-widget-content` без ошибок;
- ✅ на табе Learn присутствуют все 4 под-вкладки, переключение работает.
**Постусловия:** —

### Widgets-Spotlight-03 — Таб Notifications и «Mark all as read»
**Тип:** functional
**Предусловия:** Spotlight виден, есть непрочитанные уведомления.
**Шаги:**
1. Открыть таб **Notifications** → реакция: список уведомлений, ссылка **«Mark all as read»**, счётчик «N unread».
2. Кликнуть **«Mark all as read»** → реакция: счётчик непрочитанных обнуляется/бейдж исчезает.
**Итоговая проверка:**
- ✅ уведомления отображаются;
- ✅ после «Mark all as read» бейдж непрочитанных сбрасывается.
**Постусловия:** —
> Авто-заметка: предусловие «есть непрочитанные» нестабильно — создать уведомление через второго юзера = пошарить что-то иил написать сообещние - что проще.

---

## 8. Community — виджет ссылок

### Widgets-Community-01 — Рендер списка ссылок
**Тип:** smoke
**Предусловия:** Home page открыт, Community виден.
**Шаги:**
1. Осмотреть виджет **Community** → список ссылок (Platform Releases, Plugin releases, Macromolecules updates, Cheminformatics updates, Welcome to the Datagrok community и т.д.).
**Итоговая проверка:**
- ✅ есть минимум несколько ссылок;
- ✅ каждая ссылка ведёт на `community.datagrok.ai` (проверять `href`, а не реальный переход — внешний ресурс. но если можно проверить - будет хорошо).
**Постусловия:** —

---

## 9. Usage и Reports — виджеты UsageAnalysis (требуют прав)

> Видны только пользователям из групп **Developers / Administrators** (`canView`).
> Тестовый юзер `opavlenko+playwright` имеет роль администратора → оба виджета видны.

### Widgets-Usage-01 — Рендер Usage и «Open Usage Analysis»
**Тип:** functional
**Предусловия:** Home page под пользователем с правами; Usage виден.
**Шаги:**
1. Осмотреть виджет **Usage** → графики **Users** и **Errors** (line charts на canvas), блок **System** со ссылками (Datlas, Datlas DB, Grok Connect, Jupyter, Grok Spawner), ссылка **«Open Usage Analysis»**.
2. Кликнуть **«Open Usage Analysis»** → реакция: открывается приложение UsageAnalysis.
**Итоговая проверка:**
- ✅ присутствуют оба графика, System-ссылки и «Open Usage Analysis»;
- ✅ приложение Usage Analysis открывается без ошибок.
**Постусловия / очистка:** вернуться на Home page.

### Widgets-Reports-01 — Рендер Reports и «Open Reports»
**Тип:** functional
**Предусловия:** Home page под пользователем с правами; Reports виден.
**Шаги:**
1. Осмотреть виджет **Reports** → ссылка **«Open Reports»** и список недавних error-reports.
2. Кликнуть **«Open Reports»** → реакция: открывается вью со списком отчётов.
**Итоговая проверка:**
- ✅ присутствует ссылка «Open Reports», виджет листает недавние отчёты (или пуст);
- ✅ список отчётов открывается без ошибок.
**Постусловия / очистка:** вернуться на Home page.
> ⚠️ **Наблюдение (2026-06-08):** на dev Reports показывал 4 записи «NullError: method not found: 'clientWidth' on null» - 
> — недавние auto-error-reports от исходной ошибки (вероятно, рендер виджета в zero-width контейнере). Виджет
> работает (листает отчёты), но на саму ошибку стоит пока не обращать внимание. 

---

## 10. Permissions — гейтинг по `canView`

### Widgets-Perm-01 — Usage/Reports скрыты для пользователя без прав
**Тип:** negative
**Предусловия:** Home page под пользователем **без** ролей Developer/Administrator —
вторичный аккаунт `DATAGROK_SHARING_LOGIN` (`opavlenko+pwsharing`, storageState `.auth-sharing.json`).
**Шаги:**
1. Под non-admin аккаунтом открыть Home page.
**Итоговая проверка:**
- ✅ **Usage** и **Reports** **не отображаются**;
- ✅ видны только Spotlight и Community.
**Постусловия:** —
> Авто-заметка: убедиться, что у sharing-аккаунта нет ролей Developers/Administrators; иначе завести отдельного юзера.

### Widgets-Perm-02 — Usage/Reports видны Developer/Administrator
**Тип:** functional
**Предусловия:** Home page под основным аккаунтом `opavlenko+playwright` (имеет роль admin).
**Шаги:**
1. Открыть Home page.
**Итоговая проверка:**
- ✅ в DOM-порядке присутствуют все 4: Spotlight, Community, Usage, Reports.
**Постусловия:** —

---

## 11. Сводка по тегам (для Playwright)

- `@smoke`: Display-01, Nav-01, Spotlight-01, Community-01.
- `@functional`: Controls-01, Customize-01, Search-01, Spotlight-02/03, Usage-01, Reports-01, Perm-02.
- `@regression`: Controls-02, Customize-02.
- `@negative`: Perm-01.

**Итого: 15 кейсов.**
