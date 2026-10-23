"""Rules on tool calls: per-user paths, and the user's approval of a destructive MCP op. Every agent runs as the same
uid, so the path rules are a pattern check over the arguments: defence in depth, not isolation."""
import json
import os
import re
from typing import Any

from .protocol import JsonObject

# The /users root, not a `users` segment inside another path.
USER_DIR = re.compile(r'(?<![\w.-])/users(?![\w.-])(?:/[\w.-]+)?')
PARENT_DIR = re.compile(r'(?<![\w.])\.\.(?![\w.])')
WORKSPACE = re.compile(r'/workspace(?![\w.-])')
PROC = re.compile(r'/proc(?![\w.-])(/(?:self|thread-self|cpuinfo|meminfo|loadavg|uptime|version)(?![\w.-]))?')
STATE_DB = r'(?<![\w.-])state\.db(?![\w.-])'
# Hermes' file tools expand a leading `~` to the agent's own home, which is runtime state.
FILE_TOOLS = {'read_file', 'write_file', 'patch', 'search_files'}
EDIT_TOOLS = {'write_file', 'patch'}
# The token of a destructive MCP op's dry run, in its result text whether or not Hermes re-encoded it.
CONFIRM_TOKEN = re.compile(r'confirmToken\\?"\s*:\s*\\?"(\w+)')
CONFIRM_TOOL = 'datagrok_confirm'
NOT_APPROVED = ('The user has not approved this in the current turn. Repeat the call without `confirm` for a new dry '
                'run, ask with datagrok_confirm, and use the new token only if they approve; without that tool, it '
                'cannot be done here.')


class Approvals:
    """The MCP server's confirm token is a hash of the op's arguments that never expires, so the model could reuse
    one the user declined. A token is spendable once, and only after the user approved the confirmation card that
    followed the dry run issuing it, in the same turn; a card approves the latest dry run only."""

    def __init__(self) -> None:
        self.issued: str | None = None
        self.approved: set[str] = set()

    def check(self, tool: str, args: JsonObject) -> str | None:
        op_args = args.get('args')
        token = op_args.get('confirm') if isinstance(op_args, dict) else None
        if not tool.startswith('mcp__') or not token:
            return None
        if token not in self.approved:
            return NOT_APPROVED
        self.approved.discard(token)
        return None

    def record(self, tool: str, result: Any) -> None:
        text = result if isinstance(result, str) else json.dumps(result, default=str)
        if tool.startswith('mcp__'):
            self.issued = (CONFIRM_TOKEN.findall(text) or [self.issued])[-1]
        elif tool == CONFIRM_TOOL:
            if self.issued and re.search(r'"confirmed"\s*:\s*true', text):
                self.approved.add(self.issued)
            self.issued = None


def violation(tool: str, args: JsonObject, user_dir: str | None, homes: str, spillover: str | None = None) \
        -> str | None:
    if tool == 'skill_manage':
        return "Editing the assistant's skills is not allowed."
    text = json.dumps(args, ensure_ascii=False)
    foreign = next((m.group(0) for m in USER_DIR.finditer(text) if m.group(0) != user_dir), None)
    if foreign:
        own = f' ({user_dir})' if user_dir else ''
        return f'Access to {foreign} is blocked — you may only access your own directory{own}.'
    if PARENT_DIR.search(text):
        return 'Parent-directory (`..`) path segments are blocked — use paths inside the current directory.'
    # A user's session has no workspace but its own <user dir>/workspace.
    if user_dir and WORKSPACE.search(text.replace(f'{user_dir}/workspace', '')):
        return 'Direct access to /workspace is blocked — use cwd-relative `workspace/...` paths instead.'
    path = args.get('path')
    # Hermes leaves a large tool result in the session's own spillover directory and points the model at it.
    outside_spillover = re.sub(re.escape(spillover) + r'(?![\w.-])', '', text) if spillover else text
    if re.search(re.escape(homes) + r'(?![\w.-])|' + STATE_DB, outside_spillover) or \
            (tool in FILE_TOOLS and isinstance(path, str) and re.match(r'^~(/|$)', path)):
        return "The assistant runtime's own state is off limits."
    if any(not m.group(1) for m in PROC.finditer(text)):
        return 'Access to /proc is blocked except /proc/self.'
    return None


def absolute_edit_path(tool: str, args: JsonObject, cwd: str) -> str | None:
    """The absolute form of a relative write_file/patch path. Hermes resolves it against the session directory
    everywhere except when it prepares the edit for approval, which uses the process cwd and fails."""
    path = args.get('path')
    if tool in EDIT_TOOLS and isinstance(path, str) and path and not os.path.isabs(path):
        return os.path.normpath(os.path.join(cwd, path))
    return None
