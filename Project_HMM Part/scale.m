function VOUT = scale(VIN,v)

VOUT = [];
M = 0;
for i=1:size(VIN,1)
    if isnan(VIN(i,1))
    else
        M = M+VIN(i,1);
    end
end
M = M/size(VIN,1);

N = VIN - M(1);

big = max(abs(N));

j = 0;
while 1
    if big < 10^j
        j = j-1;
        break;
    end
    j = j + 1;
end

for i=1:size(N,1)
    VOUT = [VOUT;1/(1+exp(-v*N(i,1)/10^j))];
end

end