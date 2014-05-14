% 2-D HIO written by Po-Nan Li @ Academia Sinica 2012
function R = hio2d(Fabs, S, n, gen, rep, checker, alpha)

    % OSS module
    if numel(alpha) > 0
        oss = 1;
        x = -round((length(Fabs)-1)/2):round((length(Fabs)-1)/2);
        [X, Y] = meshgrid(x, x);
        W = exp(-0.5 .* (X./alpha).^2) .* exp(-0.5 .* (Y./alpha).^2);
        W = ifftshift(W);
    else
        oss = 0;
    end
    
    
%     cond = (rep == 1) & gen == 1;
    cond = 0;
    
    beta1 = 0.9;
    beta2 = 0.7;
    
    if numel(checker) > 0
        check = 1;
    else
        check = 0;
    end
    if sum(sum(imag(Fabs))) == 0
        ph_init = rand(size(Fabs));
%         ph_init = (ph_init + rot90(ph_init, 2) ) ./2;
        ph_init = angle(fft2(ph_init));
        F = Fabs .* exp(1j.*ph_init);
    else
        F = Fabs;
    end
    
    F0 = abs(F);

    
    previous = ifft2(F, 'symmetric');
    % ================ iterations ==================================
    for t = 1:n
        if mod(t-1, 100) == 0 && n >= 500
            disp(['step ' int2str(t)]);
        end
        rs = ifft2(F, 'symmetric'); % real space version
%         rs = rs.* S;
        cond1 = ~S | (rs<0);
        rs(cond1) = previous(cond1) - beta1 .* rs(cond1);
        previous = rs;
        if oss == 1
            rs_oss = ifft2(fft2(rs) .* W, 'symmetric');
            rs(~S) = rs_oss(~S);
        end
        F2 = fft2(rs);% .* exp(-1j.*(U+V));
        F = F0 .* exp(1j.*angle(F2));
        if check
            F(checker) = F2(checker);
        end
    end
        % ================ iterations ends here  ==================================
    R = ifft2(F, 'symmetric');
end