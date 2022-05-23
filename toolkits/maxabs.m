function a = maxabs(b)
[~,index] = max(abs(b));
a = b(index);