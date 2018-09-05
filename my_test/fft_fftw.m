
size =8;
dt=1;
df=1/(size*dt);



array_t=1:size;
array_f=fftshift(fft(ifftshift(array_t)))./(size*df);
array_f'

sum(abs(array_t).^2)*dt
2*sum(abs(array_f(4:8)).^2)*df

test=array_f(5:8);
test_2=[flip(test),0,test];
result=ifftshift(ifft(fftshift(test_2)));
plot(real(result))
hold on
plot(imag(result))
sum(abs(result).^2)*dt


a(1)=28;
a(2)=-4+1i*9.65685;
a(3)=-4+1i*4;
a(4)=-4+1i*1.65685;
a(5)=-4;


