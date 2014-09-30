#ifndef POMPYT_WRAPPER_H
#define POMPYT_WRAPPER_H

extern "C" {
    extern struct {
	double pomi[20];
	int ipom[20];
    } pompyt1_;
}
#define pompyt1 pompyt1_

extern "C" {
    extern struct {
	double pomcut[20];
    } pomcut1_;
}
#define pomcut1 pomcut1_

// POMPYT routines 

#define pomini pomini_
#define pompyt pompyt_
#define pomout pomout_

    extern "C" {
      void pomini(double*);
      void pompyt();
      void pomout();
    }

void inline call_pomini(double cmsenergy) 
{ pomini( &cmsenergy ); }

void inline call_pompyt(){ pompyt(); }
void inline call_pomout(){ pomout(); }



#endif  // POMPYT_WRAPPER_H
