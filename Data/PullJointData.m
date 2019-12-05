clear
close all

AnkleData = csvread('Ankle Angles.csv');
HipData = csvread('Hip Angles.csv');
KneeData = csvread('Knee Angles.csv');

AnkleFit = fit(AnkleData(:,1),AnkleData(:,2),'smoothingspline');
HipFit = fit(HipData(:,1),HipData(:,2),'smoothingspline');
KneeFit = fit(KneeData(:,1),KneeData(:,2),'smoothingspline');

save('C:\Users\alexj\OneDrive\Documents\UT Grad School\Brain, Body, and Robotics\BBR_Fall2019\JointData.mat')