/* Platform mentions for the compose box: users out of `grok.dapi.users`, rows and chips through the
   object handler, and the token the platform's own markup pipeline already consumes —
   `<span>#{x.<id>."<name>"}</span>` (grok_shared mixins.dart:18, rendered by the Entity markup
   handler in entity_renderers.dart). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {MessageInput, MessageInputOptions, MentionProvider} from '../../components/inputs/message-input.js';
import {dapiSource} from '../entities/dapi-source.js';
import {handlerRenderer} from '../entities/entity.js';

const SEARCH_FIELDS = ['firstName', 'lastName', 'login', 'email'];
// one renderer for the whole module: its per-collection handler cache is the point of it
const USER_RENDERER = handlerRenderer<DG.User>();

export const USER_TOKEN = '<span>#\\{x\\.([^."]+)\\."(.*?)"\\}<\\/span>';

export function userMentionProvider(): MentionProvider<DG.User> {
  return {
    trigger: '@',
    minChars: 0,
    fetch: dapiSource(() => grok.dapi.users, {filter: (q) => q ?
      SEARCH_FIELDS.map((f) => `${f} like "${q}"`).join(' or ') : undefined}),
    renderer: USER_RENDERER,
    caption: (user) => user.friendlyName ?? user.name,
    toMarkup: userMarkup,
    tokenPattern: USER_TOKEN,
    renderToken: renderUserToken,
  };
}

export function messageInput(options: Omit<MessageInputOptions, 'mentions'> = {}): MessageInput {
  return new MessageInput({...options, mentions: [userMentionProvider()]});
}

/** The real `DG.User.toMarkup()` wherever the platform supplies one; the composed fallback is the
 * same shape, down to the whitespace collapsing of grok_shared utils.dart:6-12. */
function userMarkup(user: DG.User): string {
  const own = (user as unknown as {toMarkup?: () => string}).toMarkup;
  if (typeof own === 'function')
    return own.call(user);
  const name = (user.friendlyName ?? user.name ?? '').replace(/[\n\t\r]/g, ' ');
  return `<span>#{x.${user.id}."${name}"}</span>`;
}

/** The name text stands in until the entity resolves — the placeholder-then-resolve of the
 * platform's own markup handler; a failure just keeps the name. */
function renderUserToken(match: RegExpExecArray): HTMLElement {
  const el = document.createElement('span');
  el.textContent = match[2];
  grok.dapi.users.find(match[1])
    .then((user) => {
      if (user)
        el.replaceChildren(USER_RENDERER.markup(user));
    })
    .catch(() => undefined);
  return el;
}
