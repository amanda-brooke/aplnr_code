function map = magma_colormap(n)

    if nargin < 1 || isempty(n)
        n = 256;
    end
    validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive'});
    anchors = [ 0   0   4
        12  9   38
        34  17  80
        64  15  116
        95  24  127
        123 35  130
        152 45  128
        183 55  121
        211 67  110
        235 87  96
        248 118 92
        253 154 106
        254 187 129
        253 218 156
        252 253 191] ./ 255;
    x = linspace(0, 1, size(anchors, 1));
    map = interp1(x, anchors, linspace(0, 1, double(n)), 'pchip');
    map = min(max(map, 0), 1);
end
