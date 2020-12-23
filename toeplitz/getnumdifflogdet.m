function deriv = getnumdifflogdet(Qr,dQr,N)

    circ1 = reshape(Qr,N(2),N(1));
    circ2 = reshape(dQr,N(2),N(1));
%    circ3 = reshape(d2Qr,N(2),N(1));

    circ1 = [circ1, circ1(:,(end-1):-1:2)];
    circ1 = [circ1; circ1((end-1):-1:2,:)];

    circ2 = [circ2, circ2(:,(end-1):-1:2)];
    circ2 = [circ2; circ2((end-1):-1:2,:)];

%    circ3 = [circ3, circ3(:,(end-1):-1:2)];
%    circ3 = [circ3; circ3((end-1):-1:2,:)];

    d1 = fft2(circ1); d2 = fft2(circ2); %d3 = fft2(circ3);
    dd1 = real(d2./d1);
%    dd2 = real(d3./d1);

    deriv = sum(dd1(1:N(2),1:N(1)),'all');

%    deriv2 = sum(dd2(1:N(2),1:N(1)),'all') - sum(dd1(1:N(2),1:N(1)).^2,'all');
end