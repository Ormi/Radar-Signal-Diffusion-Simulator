function make_plot(x, repetitions)
   Fs = repetitions;
   dt = 1/Fs;             

   xx = 0:dt:1;
   yy = x;

   figure

   %subplot(1,2,1);
   %plot(xx,real(yy));
   %title('Real part');
   %xlabel('Time');
   %ylabel('Power');

   %subplot(1,2,2);
   %plot(xx,imag(yy))
   plot(xx, real(yy), 'r-' , xx, imag(yy), 'b-');
   %title('Imaginary part');
   xlabel('Time');
   ylabel('Amplitude');
end
