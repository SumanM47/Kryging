function [ld, dld] = logdet(Qr,dQr,nvec)

 dim = numel(nvec);

 assert(numel(Qr) == prod(nvec(1:dim)), 'Check dimensions of Toeplitz Row');
 switch dim
     case 1
	 circ = [Qr, Qr(end:-1:2)];	
	 ff = fft(circ);
	 ff(ff<=0) = eps; 
	 dcirc = [dQr, dQr(end:-1:2)];	
	 dff = fft(dcirc);
	 dff(dff<=0) = eps;

         ld = real(sum(log(ff(1:nvec(1)))));
	 dld = real(sum(dff(1:nvec(1))./ff(1:nvec(1))));
     case 2

        circ = reshape(Qr,nvec(2),nvec(1));
        circ = [circ, circ(:,(end-1):-1:2)];
        circ = [circ; circ((end-1):-1:2,:)];
	ff = fft2(circ);
	ff(ff<=0) = eps;
        dcirc = reshape(dQr,nvec(2),nvec(1));
        dcirc = [dcirc, dcirc(:,(end-1):-1:2)];
        dcirc = [dcirc; dcirc((end-1):-1:2,:)];
	dff = fft2(dcirc);
	dff(dff<=0) = eps;

         ld = real(sum(log(ff(1:nvec(2),1:nvec(1))),'all'));
	 dld = real(sum(dff(1:nvec(2),1:nvec(1))./ff(1:nvec(2),1:nvec(1)),'all'));
	 
     otherwise
         error('Wrong dimesnion');
 end
end
