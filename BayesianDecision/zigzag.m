% zigzag.m
function output = zigzag(in)
    A = load('Zig-Zag Pattern.txt');
    output = [1:64];
    for i = 1:8
        for j = 1:8
            output(A(i, j) + 1) = in(i, j);
        end
    end
end
