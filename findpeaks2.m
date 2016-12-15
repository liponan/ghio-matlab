function [r, c, y] = findpeaks2(data)
[h w] = size(data);
r = 0;
c = 0;
y = -inf;
t = 1;
for u = 2:(h-1)
    for v = 2:(w-1)
        if (data(u,v) > data(u+1,v)) && (data(u,v) > data(u,v+1)) && (data(u,v) > data(u-1,v)) && (data(u,v) > data(u,v-1))
            r(t) = u;
            c(t) = v;
            y(t) = data(u,v);
            t = t+1;
        end
    end
end

[y, si] = sort(y,'descend');
r = r(si);
c = c(si);

end