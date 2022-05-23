function [p_coef, p_pvalues] = partialcorr_self(X,Y,Z)
%%% plot residue of X given Z v.s residue of Y give Z 
%%% plot residue of X given Y v.s residue of Z give Y
%%% plot residue of Y given X v.s residue of Y give X


vars = [X;Y;Z];
p_coef   = zeros(3,3);
p_pvalues = zeros(3,3);
for i = 1:3
    for ii = 1:3
        if i >= ii
            continue
        end
        var1 = vars(i,:);  var2 = vars(ii,:); given= vars(6-i-ii,:);
        % calculate the residue
        [rho, pval] = partialcorr(var1',var2',given');
        p_coef(i,ii)    = rho;
        p_pvalues(i,ii) = pval;
    end
end
