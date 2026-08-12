package grok_connect.table_mutation;

/**
 * The recomputed plan has destructive actions and the mutation did not set confirmDestructive —
 * maps to a structured 400 that carries the plan, so the UI renders the confirmation without a
 * second round trip (ARCHITECTURE §3.4.2).
 */
public class MutationConfirmationRequiredException extends RuntimeException {
    public final MutationPlan plan;
    public final String providerType;

    public MutationConfirmationRequiredException(MutationPlan plan, String providerType) {
        super("Destructive operation requires confirmDestructive: true");
        this.plan = plan;
        this.providerType = providerType;
    }
}
