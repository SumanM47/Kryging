function submat = toeplitzsub(ind,Qr)
    ind = cell2mat(ind);
    if ~isvector(ind)
        error('indices must be supplied in a vector format')
    end
        
    subsize = max(size(ind));
    if subsize == 1
        submat = 1;
    end
    if subsize > 1
        nall = max(size(Qr));
        submat = zeros(subsize,subsize);
        for i=1:subsize
            m = ind(i);
            if m > nall
                error('index exceeds bounds');
            end
            submat(i,i) = 1;
            bigm = find(ind > m);
            if ~isempty(bigm)
            eligind = ind(ind>m) - m;
            temp = Qr(eligind);
            submat(i,bigm) = temp(:);
            submat(bigm,i) = temp(:);
            end
        end
    end
end