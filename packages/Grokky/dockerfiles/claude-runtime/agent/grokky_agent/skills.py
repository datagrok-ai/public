"""Grokky's skills as the models see them: the index in skill_view's description, and the one redirect a turn gets
from reading Datagrok's source to the skill or tool that already covers what it is looking for."""
import json
import re
from pathlib import Path

Skill = tuple[str, str]

# The platform's source in the workspace; its documentation (workspace/help) stays open.
SOURCE = re.compile(r'(?<![\w.-])workspace/(?!help(?![\w.-]))')
SEARCH_TOOLS = {'read_file', 'search_files'}
WORD = re.compile(r'[A-Z]?[a-z]+|[A-Z]+(?![a-z])')
STOP = {'the', 'and', 'for', 'src', 'workspace', 'packages', 'api', 'get', 'set', 'with', 'from', 'datagrok'}


def index(skills_dir: str) -> list[Skill]:
    skills = []
    for path in sorted(Path(skills_dir).glob('*/SKILL.md')):
        text = path.read_text()
        head = text.split('---')[1] if text.startswith('---') else ''
        name = re.search(r'^name:\s*(.+)$', head, re.M)
        description = re.search(r'^description:\s*(.+)$', head, re.M)
        if name and description:
            skills.append((name.group(1).strip(), description.group(1).strip()))
    return skills


def _words(text: str) -> set[str]:
    return {w.lower() for w in WORD.findall(text) if len(w) > 2} - STOP


def redirect(tool: str, args: dict[str, object], skills: list[Skill]) -> str | None:
    """The message refusing a turn's first look into the platform's source, naming the skills that match it best."""
    text = json.dumps(args, ensure_ascii=False)
    if tool not in SEARCH_TOOLS or not SOURCE.search(text):
        return None
    wanted = _words(text)
    scored = sorted(((len(wanted & _words(f'{name} {desc}')), name) for name, desc in skills), reverse=True)
    names = ', '.join(name for score, name in scored[:3] if score)
    which = f'closest: {names}' if names else "the index is in skill_view's description"
    return (f'Not run: use Datagrok\'s own tooling before its source code. Load the skill that covers this with '
            f'skill_view ({which}); use the mcp__datagrok__* tools for server entities (files, functions, '
            'scripts, queries, projects, spaces, sharing) and list_view_functions for what the current view can do. '
            'If none of them covers it, repeat this call and it will run.')
