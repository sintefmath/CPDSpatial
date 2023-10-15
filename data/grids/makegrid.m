function G = makegrid(R, ncells, varargin)
    opt = merge_options(struct('angle', 6 * pi / 180, 'exponent', 1.0), varargin{:});
    
    angle = opt.angle;
    exponent = opt.exponent;
    
    spacing = linspace(0, 1, ncells+1).^exponent;
    spacing = fliplr(1-spacing)';
    spacing(1) = spacing(1) + eps; % avoid inner degenerate face
    
    nodes = spacing * R;
    %nodes = linspace(eps, R, ncells + 1)';
    
    x = cos(angle) * nodes;

    v = nodes * sin(angle);
    
    y1 = -v * sqrt(2)/2;
    y2 =  v * sqrt(2)/2;
    y3 = -v * sqrt(2)/2;
    y4 =  v * sqrt(2)/2;
    
    z1 = -v * sqrt(2)/2;
    z2 = -v * sqrt(2)/2;
    z3 =  v * sqrt(2)/2;
    z4 =  v * sqrt(2)/2;
    
    coord_x = repmat(x, 4, 1);
    coord_y = [y1;y2;y3;y4];
    coord_z = [z1;z2;z3;z4];
    
    G = cartGrid([ncells, 1, 1]);
    
    G.nodes.coords = [coord_x, coord_y, coord_z];
    
    G = computeGeometry(G);
    
end
