/* Platform user picker over the plain TypeAhead: an async source on grok.dapi.users and a
   two-line row (avatar + name + login). Nothing here is u2 machinery — swap the render and the
   same control picks projects, groups or queries. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {TypeAhead} from '../../components/inputs/typeahead.js';
import {PickInput} from '../../components/inputs/pick-input.js';
import type {InputOptions} from '../../core/input-base.js';
import {dapiSource} from '../entities/dapi-source.js';

const SEARCH_FIELDS = ['firstName', 'lastName', 'login', 'email'];
const AVATAR_HUES = 6;

/** Smart-filter `like` wraps the value in `%…%` server-side, so the sanitized substring goes
 * in as is. */
export function userInput(placeholder = 'User…'): TypeAhead<DG.User> {
  return new TypeAhead<DG.User>({
    source: dapiSource(() => grok.dapi.users,
      {filter: (q) => q ? SEARCH_FIELDS.map((f) => `${f} like "${q}"`).join(' or ') : undefined}),
    itemText: (user) => user.friendlyName,
    keyOf: (user) => user.id,
    render: renderUser,
    placeholder,
    openOnFocus: true,
    maxHeight: 320,
  });
}

/** A `User` reference's editor: the value is the user id, the box the user picker. */
export function userPick(options: InputOptions<string | null> = {}): PickInput<DG.User> {
  return new PickInput<DG.User>({...options, typeAhead: () => userInput(''), idOf: (user) => user.id,
    resolve: (id) => grok.dapi.users.find(id)});
}

function renderUser(user: DG.User): HTMLElement {
  const row = document.createElement('div');
  row.className = 'u2-typeahead-user';
  const text = document.createElement('div');
  text.className = 'u2-typeahead-user-text';
  const name = line('u2-typeahead-user-name', user.friendlyName);
  name.dataset.u2Part = 'item-label';
  text.append(name);
  const secondary = user.login || user.email;
  if (secondary && secondary !== user.friendlyName)
    text.append(line('u2-typeahead-user-secondary', secondary));
  row.append(avatar(user), text);
  return row;
}

/** `User.picture` is typed `string | object` but the Dart binding renders the avatar element
 * (grok_api.dart: `User_Get_Picture` → `UserMeta.renderIcon`); initials cover anything else. */
function avatar(user: DG.User): HTMLElement {
  const picture = user.picture;
  if (picture instanceof HTMLElement) {
    const host = line('u2-avatar u2-avatar-picture');
    host.append(picture);
    return host;
  }
  return line(`u2-avatar u2-avatar-${hue(user.login)}`, initials(user.friendlyName));
}

function initials(name: string): string {
  const local = name.includes('@') ? name.substring(0, name.indexOf('@')) : name;
  const capitals = local.match(/[A-Z]/g);
  return capitals ? capitals.slice(0, 2).join('') : local.substring(0, 2).toUpperCase();
}

function hue(login: string): number {
  let hash = 0;
  for (let i = 0; i < login.length; i++)
    hash = (hash * 31 + login.charCodeAt(i)) % 9973;
  return hash % AVATAR_HUES + 1;
}

function line(cls: string, text?: string): HTMLElement {
  const el = document.createElement('div');
  el.className = cls;
  if (text !== undefined)
    el.textContent = text;
  return el;
}
