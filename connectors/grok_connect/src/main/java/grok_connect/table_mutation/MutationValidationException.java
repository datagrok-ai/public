package grok_connect.table_mutation;

import java.util.List;

/** Invalid mutation payload (missing lists, empty WHERE without allowFullTable, unknown dg type) — maps to a structured 4xx. */
public class MutationValidationException extends RuntimeException {
    /** Hard-refusal details ({path, code, message, count}) for live-data pre-check failures; null otherwise. */
    public List<DestructiveAction> refusals;

    public MutationValidationException(String message) {
        super(message);
    }

    public MutationValidationException(String message, List<DestructiveAction> refusals) {
        super(message);
        this.refusals = refusals;
    }
}
