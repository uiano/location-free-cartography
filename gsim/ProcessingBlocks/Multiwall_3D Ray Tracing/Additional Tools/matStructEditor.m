clear
clc

[structFileName,structFileAddress] = uigetfile('.mat');
load([structFileAddress,structFileName]);