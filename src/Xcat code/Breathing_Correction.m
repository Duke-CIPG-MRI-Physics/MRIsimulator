close all;clear;clc;

AP_curve = readmatrix('C:\Users\Roberto\Documents\XCAT\XCAT_latest\ap_curve.dat');

t = AP_curve(:,1);
AP_curve = AP_curve(:,2);

AP_curve_fit = fit(t,AP_curve,'sin3');

hold on
plot(t,AP_curve)
plot(f)

