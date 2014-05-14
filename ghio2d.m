% 2-D Guided-HIO written by Po-Nan Li @ Academia Sinica 2012
% Reference: Chen et al., Phys. Rev. B. 76, 064113 (2007).
function R = ghio2d(Fabs, S, n, gen, rep, checker, alpha)
    R = zeros(size(Fabs,1), size(Fabs,2), gen);
    G = zeros(size(Fabs,1), size(Fabs,2), rep);
    efs = zeros(rep, 1);
    sh = floor(sqrt(numel(find(S==1)))/2);
    for g = 1:gen
        disp(['in generation #' int2str(g)]);
        if g == 1
            parfor r = 1:rep
                G(:,:,r) = hio2d(Fabs, S, n, checker, alpha);
            end
        else
            parfor r = 1:rep
                G(:,:,r) = hio2d(fft2(G(:,:,r)), S, n, checker, alpha);
            end
        end
        FG = fft2(G);
        for r = 1:rep
            efs(r) = ef(Fabs, FG(:,:,r), checker);
        end
        [min_ef, min_ind] = min(efs);
        disp(['replica #' int2str(min_ind) ' selected (EF = ' num2str(min_ef) ')']);
        GM = G(:,:,min_ind);
        GM(GM<0) = 0;
        parfor r = 1:rep
            G(:,:,r) = align2(GM, G(:,:,r), sh);
            G(:,:,r) = sqrt( G(:,:,r) .* GM );
        end
        R(:,:,g) = GM;
    end
   


end