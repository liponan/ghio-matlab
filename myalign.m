function B = myalign(A, B)
    xcorr = fftshift( ifft2( fft2(A) .* conj(fft2(B)) ) );
    [cr1, cc1, y1] = findpeaks2(xcorr);
    xcorr = fftshift( ifft2( fft2(A) .* conj(fft2( rot90(B,2) )) ) );
    [cr2, cc2, y2] = findpeaks2(xcorr);
    if y2(1) > y1(1)
        B = rot90(B, 2);
        cr2 = cr2(1) - round(size(A,1)/2);
        cc2 = cc2(1) - round(size(A,2)/2);
        B = circshift(B, [cr2, cc2]);
        disp(['myalign: inverted and shifted (' int2str(cr2) ', ' int2str(cc2) ').']);
    else
        cr1 = cr1(1)- round(size(A,1)/2);
        cc1 = cc1(1) - round(size(A,2)/2);
        B = circshift(B, [cr1, cc1]);
        disp(['myalign: shifted (' int2str(cr1) ', ' int2str(cc1) ').']);
    end           
    
end