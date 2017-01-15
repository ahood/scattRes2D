function show_error(obj,r)
    E = obj.compute_error(r);
    contour(r.X,r.Y,log10(E)); colorbar;
end