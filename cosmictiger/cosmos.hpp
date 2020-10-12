#pragma once

struct cosmos {
	double t, tau, a, adot;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & t;
		arc & tau;
		arc & a;
		arc & adot;
	}
	cosmos();
	void advance_to_time(double, int = 100);
	void advance_to_scalefactor(double);
};
