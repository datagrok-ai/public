---
title: "Managed settings"
mdx:
  format: mdx
sidebar_position: 4
keywords:
  - settings
  - defaults
  - lock
  - group settings
  - settings override
---

Administrators and group admins can define default client settings for a
[group](users-and-groups.md#groups) or for everyone on the platform, and
optionally **lock** individual settings or whole settings pages so that group
members can't change them. Typical uses: enforcing a corporate color palette,
turning off beta features for everyone, or pre-configuring formatting defaults
for a team.

Shared defaults are managed from the same place as personal settings: on the
**Sidebar**, click **Settings** (<FAIcon icon="fa-solid fa-gear"/>).

## Who can manage shared settings

| You have...                                     | You can edit...                                        |
|-------------------------------------------------|--------------------------------------------------------|
| The `EditSettings` [global permission](access-control.md) | The **All users** scope and any group's defaults |
| Admin membership in a group                     | That group's defaults (and its subgroups')             |
| Neither                                         | Only your personal settings (no scope toggle appears)  |

## Set defaults for a group

1. On the **Sidebar**, click **Settings** (<FAIcon icon="fa-solid fa-gear"/>).
1. In the ribbon, use the scope toggle: the **laptop** icon shows your personal
   settings, the **group** icon switches to shared defaults.
1. Pick a scope from the dropdown next to the toggle: **All users**, or search
   for any group or role you're allowed to edit.
1. Edit values on any client settings page (General, Windows, Formats, Colors,
   and so on).
1. Click **Save and apply**.

Edits are staged per scope, so you can move between pages and scopes before
saving. Only the settings you actually change become group defaults — everything
else stays untouched.

Members pick up the new defaults the next time they sign in or refresh the
platform.

## Lock settings

By default, a group default is just that — a default. A member who changes the
setting keeps their own value. To prevent that, lock the setting:

* **Lock one setting**: In a group or All-users scope, each setting has a lock
  toggle next to it. Click it to lock or unlock the setting, then click
  **Save and apply**.
* **Lock a whole page**: Click the page-lock button in the ribbon to lock every
  setting on the open page. A per-setting lock set explicitly always wins over
  the page lock, so you can lock a page and still leave individual settings
  open, or the other way around.
* **Lock without a value**: You can lock a setting without giving it a value.
  Members keep whatever value they currently have, but can no longer change it.

For members, a locked setting appears disabled with a lock icon. Hovering over
it shows who manages it (**Managed by your administrator**).

### How values and locks are resolved

For each setting, a member's effective value is, in order:

1. Their personal value — unless the setting is locked for them.
1. The default from the nearest group they belong to that defines one (a
   direct group beats a more distant parent group).
1. The **All users** default.
1. The built-in platform default.

Locks follow the same "nearest wins" rule: a group-level lock or unlock
overrides an **All users** lock for that group's members.

:::note

Locks control the platform user interface and the defaults applied at sign-in.
They are a management convenience, not a security boundary — don't use them to
protect sensitive values.

:::

## Restore defaults

In your personal scope:

* Every setting you've changed shows a **Restore default** control next to it.
  It restores the value your administrators defined for you or, if there is
  none, the built-in default.
* The **Reset this section to defaults** ribbon button resets every unlocked
  setting on the open page.

## Platform-managed settings pages

On instances managed by an external administration console (such as Datagrok
SaaS tenants), the server-side **Admin** and **Users** settings pages are
read-only and show a **Managed by your administrator** banner with an
**Open admin console** link. Change those settings in the admin console —
they're applied to the instance from there.

See also:

* [Users and groups](users-and-groups.md)
* [Access control](access-control.md)
