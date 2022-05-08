package utilities.extetions

inline fun <reified T> T?.orIfNull(input: () -> T): T {
    return this ?: input()
}