function Bd = dilation(B,SE)

Bd = logical(convn(B,SE,'same'));