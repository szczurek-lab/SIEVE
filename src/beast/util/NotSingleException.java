package beast.util;

import beast.core.Description;

@Description("Throw this exception when exactly one object is required but more or less found.")
public class NotSingleException extends Exception {

    private static final long serialVersionUID = -5264539886225844179L;

    public NotSingleException() {
    }

    public NotSingleException(String message) {
        super(message);
    }

    public NotSingleException(String message, Throwable cause) {
        super(message, cause);
    }

    public NotSingleException(Throwable cause) {
        super(cause);
    }

}
