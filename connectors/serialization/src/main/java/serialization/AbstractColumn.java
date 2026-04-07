package serialization;

import java.util.HashMap;
import java.util.Map;

public abstract class AbstractColumn<T> implements Column<T> {
    public static final int DEFAULT_INIT_SIZE = 100;
    private final Map<String, String> tags = new HashMap<>();
    protected String name;
    protected int length = 0;
    protected int initColumnSize;

    public AbstractColumn(String name) {
        this(name, DEFAULT_INIT_SIZE);
    }

    public AbstractColumn(String name, int initColumnSize) {
        this.name = name;
        this.initColumnSize = initColumnSize;
    }

    @Override
    public String getName() { return name; }

    @Override
    public void setName(String name) { this.name = name; }

    @Override
    public int getLength() { return length; }

    @Override
    public int getInitColumnSize() { return initColumnSize; }

    @Override
    public void setInitColumnSize(int size) { this.initColumnSize = size; }

    @Override
    public Map<String, String> getTags() { return tags; }
}
