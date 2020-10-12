#pragma once



struct cosmos {
	double t, tau, a, adot;
	cosmos();
	void advance_to_time(double,int=100);
	void advance_to_scalefactor(double);
};
