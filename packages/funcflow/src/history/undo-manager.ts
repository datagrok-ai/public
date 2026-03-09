/** Command Pattern undo/redo manager */

export interface Command {
  execute(): void;
  undo(): void;
  description: string;
}

export class UndoManager {
  private undoStack: Command[] = [];
  private redoStack: Command[] = [];
  private maxHistory = 100;

  /** Execute a command and add it to the undo stack */
  execute(command: Command): void {
    command.execute();
    this.undoStack.push(command);
    this.redoStack = [];
    if (this.undoStack.length > this.maxHistory)
      this.undoStack.shift();
  }

  /** Record a command that was already executed (for external events) */
  record(command: Command): void {
    this.undoStack.push(command);
    this.redoStack = [];
    if (this.undoStack.length > this.maxHistory)
      this.undoStack.shift();
  }

  /** Undo the last command */
  undo(): boolean {
    const command = this.undoStack.pop();
    if (!command) return false;
    command.undo();
    this.redoStack.push(command);
    return true;
  }

  /** Redo the last undone command */
  redo(): boolean {
    const command = this.redoStack.pop();
    if (!command) return false;
    command.execute();
    this.undoStack.push(command);
    return true;
  }

  /** Check if undo is available */
  canUndo(): boolean {return this.undoStack.length > 0;}

  /** Check if redo is available */
  canRedo(): boolean {return this.redoStack.length > 0;}

  /** Clear all history */
  clear(): void {
    this.undoStack = [];
    this.redoStack = [];
  }
}
