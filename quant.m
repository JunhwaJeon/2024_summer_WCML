
x=-1:0.01:1;
a=3;
b=round((2^(a-1))*x)/(2^(a-1));

plot(x,b,'x');
title('Quantization function (alpha=3)');
xlabel('Weight [-1,1]');
ylabel('Qauntized weight');
grid on;