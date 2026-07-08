---
feature: chat
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: []
realized_as:
  - projects-chat-collaboration-spec.ts
related_bugs: []
---

# Comments / Chat collaboration on Projects

Verifies that two collaborators can use the Chat (Comments) panel on a shared project to post and reply to comments, and that the other user sees new comments, notifications, and mentions correctly. Also checks that comment history persists across logins and that users can only edit or delete their own comments.

**Preconditions:**
- Two users available: **User A** (author) and **User B** (collaborator), each logged in to a separate browser session.
- A saved project exists that both users can access (User A shares it during the test).

1. Open chat panel:
- As **User A**, open **Browse > Projects** and select a saved project.
- In the Context Panel expand the **Chat** (Comments) section.
- Expected Result: The chat area with a message input field is shown; empty if there are no comments yet.

2. Post a comment:
- As **User A**, type `First comment from A` and send it (Enter or Send button).
- Expected Result: The comment appears in the feed with author (User A), avatar and timestamp.

3. Share the project:
- As **User A**, click **Share** and grant **User B** View/Edit access.
- Expected Result: User B gets access; a sharing notification is generated.

4. Collaborator sees the comment:
- As **User B**, open the same project and expand the **Chat** section.
- Expected Result: User B sees User A's comment with correct author and time.

5. Notification:
- As **User B**, check the notifications / inbox for a new comment or shared-project alert.
- Expected Result: A notification is present and navigates to the project/comment.

6. Reply:
- As **User B**, post `Reply from B`.
- Expected Result: The reply is added to the feed in chronological order after User A's comment.

7. Author sees the reply:
- As **User A**, refresh / reopen the project chat panel.
- Expected Result: User A sees User B's reply; messages are chronological and authors are correct.

8. Mention (if supported):
- As **User A**, post a comment with a mention, e.g. `@UserB please review`.
- Expected Result: The mention is highlighted; User B receives a personal mention notification.

9. Edit / delete own comment (if available):
- As **User A**, edit or delete one of your own comments.
- Expected Result: The comment is updated/removed; the change is reflected for User B after refresh.

10. Permissions on others' comments:
- As **User B**, try to delete User A's comment.
- Expected Result: The action is unavailable / forbidden — a user can only manage their own comments.

11. Persistence:
- Re-login as both users and reopen the project.
- Expected Result: The full chat history is preserved on the server; messages and their order are not lost.

**Postconditions / Cleanup:**
- Delete the test comments.
- Revoke User B's access to the project if it was shared only for this test.

---
{
  "order": 1
}
