function B = myalign(A, B)
    xcorr = fftshift( ifft2( fft2(A) .* conj(fft2(B)) ) );
    [cr1, cc1, y1] = findpeaks2(xcorr);
    xcorr = fftshift( ifft2( fft2(A) .* conj(fft2( rot90(B,2) )) ) );
    [cr2, cc2, y2] = findpeaks2(xcorr);
    if y2 > y1
        B = rot90(B, 2);
        cr2 = cr2 - round(size(A,1)/2);
        cc2 = cc2 - round(size(A,1)/2);
        B = circshift(B, [-cr2, -cc2]);
    else
        cr1 = cr1 - round(size(A,1)/2);
        cc1 = cc1 - round(size(A,1)/2);
        B = circshift(B, [-cr1, -cc1]);
    end           
    
end