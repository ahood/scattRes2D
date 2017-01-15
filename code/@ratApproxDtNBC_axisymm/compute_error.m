function E = compute_error(obj,r)
    % given a rectangle object, look at error on it
    E = 0*r.Z - Inf; % set up error on grid
%     for jj = 1:obj.dtn.Ns + 1 % look at nonpositive ones (all that's needed)
    for jj = 1:obj.dtn.Nt
        ratfn = obj.ratf(jj);
        E = max(E,ratfn.compute_error(r));
    end 
end