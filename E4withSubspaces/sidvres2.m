function f = sidvres2(Beta, Y, X)

Beta = Beta.^2;
f = sum((Y-X*Beta).^2);

