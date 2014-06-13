% 2-D Guided-HIO written by Po-Nan Li @ Academia Sinica 2012
% Reference: Chen et al., Phys. Rev. B. 76, 064113 (2007).
function [R, G, efs] = ghio2d(Fabs, S, n, gen, rep, checker, alpha)
    R = zeros(size(Fabs,1), size(Fabs,2), gen);
    G = zeros(size(Fabs,1), size(Fabs,2), rep);
    FG = zeros(size(Fabs,1), size(Fabs,2), rep);
    efs = zeros(gen, rep);
    sh = floor(sqrt(numel(find(S==1)))/2); % maximun shift
    GM = zeros(size(Fabs));

    for g = 1:gen
        disp(['in generation #' int2str(g)]);
        if g == 1
            parfor r = 1:rep
                G(:,:,r) = hio2d(Fabs, S, n, checker, alpha);
                FG(:,:,r) = fft2( G(:,:,r) );
                efs(g,r) = ef(Fabs, FG(:,:,r), checker);
            end
        else
            parfor r = 1:rep
                % make new template
                G(:,:,r) = myalign( GM, G(:,:,r) );
                G(:,:,r) = sign( GM ) .* sqrt( abs(G(:,:,r) .* GM) );
                % run independent HIO
                G(:,:,r) = hio2d(fft2(G(:,:,r)), S, n, checker, alpha);
                % additional treat to FG
                FG(:,:,r) = G(:,:,r);
                FG(:,:,r) = FG(:,:,r) .* S;
                FG(:,:,r) = fft2( FG(:,:,r) );
                efs(g,r) = ef(Fabs, FG(:,:,r), checker);
            end
        end
        % find the best replica
        [min_ef, min_ind] = min(efs(g,:));
        disp(['replica #' int2str(min_ind) ' selected (EF = ' num2str(min_ef) ')']);
        GM = G(:,:,min_ind);
        R(:,:,g) = GM;
        disp('===========================================================');
    end
   


end