function err = er(A,R,S)
    if numel(S) > 0
        err = sum(sum( abs( abs(A(S)) - abs(R(S)) ) ));
        err = err ./ sum(sum( abs(A(S)) ));
    else
        err = sum(sum( abs( abs(A) - abs(R) ) ));
        err = err ./ sum(sum( abs(A) ));
    end
        
end