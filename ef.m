function err = ef(Fabs, F, checker)
    err = sum(sum( abs( abs(Fabs(~checker)) - abs(F(~checker)) ) ));
    err = err ./ sum(sum( abs(Fabs(~checker)) ));
end