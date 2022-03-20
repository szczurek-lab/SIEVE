package beast.util;

public class NoFileHeaderFoundException extends Exception {

    private static final long serialVersionUID = -4302008734551291520L;

    public NoFileHeaderFoundException() {
    }

    public NoFileHeaderFoundException(String message) {
        super(message);
    }

    public NoFileHeaderFoundException(String message, Throwable cause) {
        super(message, cause);
    }

    public NoFileHeaderFoundException(Throwable cause) {
        super(cause);
    }

}
