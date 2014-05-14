function R = hio3d(Fabs, S, n, checker)


    
    beta1 = 0.9;
    
    if numel(checker) > 0
        check = 1;
    else
        check = 0;
    end

    if sum(sum(sum(abs(imag(Fabs))))) == 0
        ph_init = rand(size(Fabs));
    %         ph_init = (ph_init + rot90(ph_init, 2) ) ./2;
        ph_init = angle(fftn(ph_init));
        F = Fabs .* exp(1j.*ph_init);
    else
        F = Fabs;
    end

    
    F0 = abs(F);

    
    previous = ifftn(F, 'symmetric');
    % ================ iterations ==================================
    for t = 1:n

            disp(['step ' int2str(t)]);

        rs = ifftn(F, 'symmetric'); % real space version
%         rs = rs.* S;
        cond1 = ~S | (rs<0);
        rs(cond1) = previous(cond1) - beta1 .* rs(cond1);
        previous = rs;

        F2 = fftn(rs);% .* exp(-1j.*(U+V));
        F = F0 .* exp(1j.*angle(F2));
        if check
            F(checker) = F2(checker);
        end
    end
        % ================ iterations ends here  ==================================
    R = ifftn(F, 'symmetric');
end