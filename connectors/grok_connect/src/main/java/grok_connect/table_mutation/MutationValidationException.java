package grok_connect.table_mutation;

/** Invalid mutation payload (missing lists, empty WHERE without allowFullTable, unknown dg type) — maps to a structured 4xx. */
public class MutationValidationException extends RuntimeException {
    public MutationValidationException(String message) {
        super(message);
    }
}
