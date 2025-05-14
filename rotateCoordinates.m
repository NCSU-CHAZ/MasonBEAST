function[RX, RY] = rotateCoordinates(X, Y, LocalX, LocalY, theta)
    % ROTATECOORDINATES Rotates the coordinates (X, Y) around the origin
    % defined by (LocalX, LocalY) by an angle theta.
    %
    % Inputs:
    %   - X      : nx1 or 1xn vector of X coordinates
    %   - Y      : nx1 or 1xn vector of Y coordinates
    %   - LocalX : scalar value for the X coordinate of the local origin
    %   - LocalY : scalar value for the Y coordinate of the local origin
    %   - theta  : scalar value of the rotation angle in degrees
    %
    % Outputs:
    %   - Xout   : nx1 or 1xn vector of rotated X coordinates
    %   - Yout   : nx1 or 1xn vector of rotated Y coordinates
    
    % Translate coordinates to the origin
    Xo = X - LocalX;
    Yo = Y - LocalY;
    
    % Perform the rotation
    RY = Yo .* cosd(theta) + Xo .* sind(theta);
    RX = Xo .* cosd(theta) - Yo .* sind(theta);

    % Translate to local coordinates
    RX = RX-21;