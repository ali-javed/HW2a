%Script file to call comparedervis.m function and initiate the deravitive
%function variables

%Created by: Zoe Koch, Ali Javed, Feb-13-2018
%Class Project1 CS302 - Modeling Complex Systems
%%
%house keeping
clear all
close all
clc


%% declare second and third order deravitive to pass to comparedervis
negSin = @(x) - sin(x)
negCos = @(x) - cos(x)
%call comparedervis
comparederivs([pi/4],@sin, @cos, negSin, negCos)

