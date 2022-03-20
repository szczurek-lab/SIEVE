package beast.util;

public class NoMatchLengthException extends Exception {

    private static final long serialVersionUID = -2395596453505590339L;

    public NoMatchLengthException() {
    }

    public NoMatchLengthException(String message) {
        super(message);
    }

    public NoMatchLengthException(String message, Throwable cause) {
        super(message, cause);
    }

    public NoMatchLengthException(Throwable cause) {
        super(cause);
    }

}
