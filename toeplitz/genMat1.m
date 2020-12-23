function out = genMat1(xmin,xmax,nvec,nu,range,rval,theta)
	if nu <= 0 | range <= 0 | rval <= 0 | rval > 1
		error('improper input');
	end

	if nu == 0.5
		k = @(r) expcov(r,range,rval);
	else
		k = @(r) matern(r,nu,range,rval);
	end


	nvec0 = nvec;
 	xmax0 = xmax;
	posdef = false;
	maxcount = 5;
	counter = 0;

	while (~posdef & (counter < maxcount))
		Qr = createrow(xmin,xmax0,nvec0,k,theta);
	        circ = reshape(Qr,nvec0(2),nvec0(1));
	        circ = [circ,circ(:,(end-1):-1:2)];
	        circ = [circ;circ((end-1):-1:2,:)];
	        [n1,n2] = size(circ);
		ff = real(fft2(circ));
		if (min(min(ff))<=0)
			counter = counter + 1;
			xmax0 = xmin + ((xmax0 - xmin).*(2.*nvec0 - 1)/(nvec0 - 1));
			nvec0 = 2.*nvec0;
		else
			posdef = true;
		end
	end
	if counter == maxcount
		error('Could not find positive definite embedding for Q');
        end

	F = ifft2(sqrt(ff).*fft2(randn([n1,n2])));
	out = real(F(1:nvec(2),1:nvec(1)));
end