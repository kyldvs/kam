package com.kyldvs.kam.matrix.misc;

public class InvalidDimensionsException extends RuntimeException {

	private static final long serialVersionUID = 1L;
	
    public InvalidDimensionsException() {
        super();
    }

    public InvalidDimensionsException(String message) {
        super(message);
    }

    public InvalidDimensionsException(String message, Throwable cause) {
        super(message, cause);
    }

    public InvalidDimensionsException(Throwable cause) {
        super(cause);
    }
}
