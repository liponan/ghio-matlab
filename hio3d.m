% 3-D HIO written by Po-Nan Li @ Academia Sinica 2012
function R = hio3d(Fabs, S, n, varargin)

    beta1 = 0.9; % feedback parameter
    
    if isempty(varargin)
        unknown = false(size(Fabs));
    else
        unknown = varargin{1};
    end

    % generate initial random phases
    if sum(abs(imag( Fabs(:) ))) == 0
        ph_init = rand(size(Fabs));
        ph_init = angle(fftn(ph_init));
        F = Fabs .* exp(1j*ph_init);
    else
        F = Fabs;
    end

    F0 = abs(F);
    previous = ifftn(F, 'symmetric');
    % ================ iterations ==================================
    for t = 1:n
        if mod(t-1, 100) && n > 499
            disp(['step ' int2str(t)]);
        end

        rs = ifftn(F, 'symmetric'); % real space version
        % apply the support
        cond1 = ~S | (rs<0);
        rs(cond1) = previous(cond1) - beta1 .* rs(cond1);
        previous = rs; % update the buffer

        F2 = fftn(rs);
        F = F0 .* exp(1j*angle(F2));
        F(unknown) = F2(unknown);
    end
        % ================ iterations ends here  ==================================
    R = ifftn(F, 'symmetric');
end