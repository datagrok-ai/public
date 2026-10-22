"""python -m grokky_agent: one chat session. The first stdin line is the init command (credentials and settings, so no
secret is in this process's environment); every later line is a command, every stdout line an event. With --check it
instead makes one short call to the model, prints a `checked` event and exits: the runtime's connection test."""
import ctypes
import json
import os
import queue
import sys
import threading

PR_SET_DUMPABLE = 4
# Before anything else: the other agents run as the same uid and must not read this process's memory or fds.
ctypes.CDLL(None).prctl(PR_SET_DUMPABLE, 0, 0, 0, 0)

from pathlib import Path  # noqa: E402

from . import config, hermes, providers, session  # noqa: E402
from .protocol import Command, ContinueCommand, Event, InitCommand, TurnCommand  # noqa: E402

TurnQueue = queue.Queue[tuple[TurnCommand | ContinueCommand, session.Turn]]


def run_turns(s: session.Session, turns: TurnQueue) -> None:
    while True:
        cmd, turn = turns.get()
        try:
            if cmd['type'] == 'turn':
                s.start_turn(cmd, turn)
            else:
                s.continue_turn(cmd, turn)
        except Exception as e:
            s.emit({'type': 'error', 'message': f'{type(e).__name__}: {e}'})


def main() -> int:
    out = os.fdopen(os.dup(1), 'w', buffering=1)
    os.dup2(2, 1)  # anything Hermes prints goes to stderr; stdout carries only events
    lock = threading.Lock()

    def emit(event: Event) -> None:
        with lock:
            out.write(json.dumps(event, default=str) + '\n')

    init: InitCommand = json.loads(sys.stdin.readline())
    try:
        provider = providers.from_credentials(init['credentials'])
    except providers.ProviderError as e:
        emit({'type': 'fatal', 'message': str(e)})
        return 1
    settings = init['settings']
    config.write_home(Path(os.environ['HERMES_HOME']), settings)
    session.current = s = session.Session(emit, provider, settings['homes'], hermes.open_store(),
                                          hermes.tool_registry(), hermes.spillover_dir())
    hermes.load_plugins()
    if not s.hooks_registered:
        emit({'type': 'fatal', 'message': 'the Grokky plugin did not load'})
        return 1
    hermes.send_keys_by_type()
    if '--check' in sys.argv:
        emit({'type': 'checked', 'error': s.check()})
        return 0
    if settings['warm']:
        s.build_agent()
    emit({'type': 'ready'})

    # Turns run on one worker: a cancel and the browser's answers must reach a turn while it runs.
    turns: TurnQueue = queue.Queue()
    threading.Thread(target=run_turns, args=(s, turns), daemon=True).start()
    for line in sys.stdin:
        cmd: Command = json.loads(line)
        if cmd['type'] == 'cancel':
            s.cancel()
        elif cmd['type'] == 'browser_result':
            s.browser_answer(cmd)
        else:
            turns.put((cmd, s.receive(cmd)))
    return 0


sys.exit(main())
