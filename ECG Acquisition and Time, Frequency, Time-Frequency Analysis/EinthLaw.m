function [Lead2_estimated, abs_diff, mse] = EinthLaw(data1,data2,data3)
Lead2_estimated=data1 + data3;
B= (data2-Lead2_estimated);
abs_diff= abs(B);
mse=(rms(abs_diff)).^2;

end

