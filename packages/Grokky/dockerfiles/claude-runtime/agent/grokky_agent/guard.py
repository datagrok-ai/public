"""Per-user path rules on tool arguments. Every agent runs as the same uid, so this is a pattern check over the
arguments: defence in depth, not isolation."""
import json
import os
import re

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
