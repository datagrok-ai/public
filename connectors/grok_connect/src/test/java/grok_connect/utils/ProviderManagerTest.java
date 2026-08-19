package grok_connect.utils;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.NullSource;
import org.junit.jupiter.params.provider.ValueSource;
import java.util.Arrays;
import java.util.LinkedHashSet;

/**
 * providers.conf is written by the Dockerfile with a trailing newline (printf '%s\n'), so the
 * main image's default empty allowlist arrives as "\n" — it must parse as "all providers"
 * (null), not an empty set that de-advertises everything (the 2.8.3 regression).
 */
class ProviderManagerTest {

    @DisplayName("All-providers shapes: null, empty, whitespace-only, wildcard")
    @ParameterizedTest
    @NullSource
    @ValueSource(strings = {"", "\n", " \n", "*", " * \n"})
    public void allProviders(String value) {
        Assertions.assertNull(ProviderManager.parseAllowlist(value));
    }

    @Test
    public void allowlistTrimsTokensAndTrailingNewline() {
        Assertions.assertEquals(new LinkedHashSet<>(Arrays.asList("Neptune", "Impala")),
                ProviderManager.parseAllowlist("Neptune, Impala\n"));
    }
}
