#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <time.h>	
#include <algorithm>	//std::copy
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

struct Walker_info{
	double x;
	double weight;
	int number;
	double Qavg;
};

double fRand(double fMin, double fMax)
{
	double f = (double)rand()/RAND_MAX;
	return fMin + f*(fMax-fMin);
}

// this function takes the initial value x and modify it
// the OU equation is solved by modified Runge-Kutta scheme
// and returns the cumulative current Q 
double propagation(double& x, double tint, gsl_rng* rng)
{
        double h = 0.001;
        double k1,k2;
        double xnew = 0.0;
        double gamma = 0.1;
        double epsilon = 1.0;
        int N = tint/h;
        double sum = 0.0;
        double a;
        int i;
        for(i=0;i<N;i++){
                a = gsl_ran_gaussian(rng,1);
                k1 = (-1)*gamma*x;
                k2 = (-1)*gamma*(x+h*k1+sqrt(epsilon*h)*a);
                xnew = x + 0.5*h*(k1+k2)+sqrt(epsilon*h)*a;
                sum += xnew;
                x = xnew;
        }
        sum = sum*h;
        return sum;
}

int main()
{
        gsl_rng* rng;
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,time(NULL));

	srand(time(0));

	int Nw = 2000;
	int multicount;
	double beta = 0.05;
	double tint = 0.1, tobs = 1000.0;
	int Ntint = tobs/tint;
	double weightsum = 0.0;
	int walkersum = 0;
	int i,j,k,s;
	double q,Q,Q2,sigma,phi;
	vector<double> mean,var,ldf,multiplicity;

	//each walker is associated with its position, weight, number, Qavg
	vector<Walker_info> walker,dummy;
	Walker_info initialize = {0.0, 0.0, 0, 0.0};

	for(i=0;i<Nw;i++){
		walker.push_back(initialize);
	}

	for(i=0;i<Ntint;i++){
		walkersum = 0;
		weightsum = 0.0;
		multicount = Nw;
		for(j=0;j<Nw;j++){
			q = propagation(walker[j].x,tint,rng);
			walker[j].Qavg += q;
			walker[j].weight = exp(-1*beta*q);
			weightsum += walker[j].weight;
		}
		for(j=0;j<Nw;j++){
			walker[j].number = floor(Nw*walker[j].weight/weightsum + fRand(0,1));
			walkersum += walker[j].number;
		}
		if(walkersum < Nw){
			while(walkersum < Nw){
				walker[fRand(0,Nw)].number += 1;
				walkersum += 1;
			}	
		}
		if(walkersum > Nw){
			while(walkersum > Nw){
				s = floor(distributionNw(generator));
				printf("%d\n",s);
				if(walker[s].number>0){
					walker[s].number -= 1;
					walkersum -= 1;
				}
			}
		}		 	
		for(j=0;j<Nw;j++){
			if(walker[j].number==0){
				multicount--;
			}
			else{
				for(k=0;k<walker[j].number;k++){
					dummy.push_back(walker[j]);
				}
			}
		}
		if(dummy.size()==Nw){
			walker.erase(walker.begin(),walker.end());
			copy(dummy.begin(),dummy.end(),walker.begin());
			dummy.erase(dummy.begin(),dummy.end());	
		}

		//evaluating average observable and large deviation function
		for(j=0;j<Nw;j++){
			Q += walker[j].Qavg;
			Q2 += walker[j].Qavg * walker[j].Qavg;
			phi += exp(-1*beta*walker[j].Qavg);
		}
		Q = Q/Nw;
		Q2 = Q2/Nw;
		sigma = (Q2-Q*Q)/tobs;
		Q = Q/tobs;
		phi = log(phi/Nw)/tobs;
		mean.push_back(Q);
		var.push_back(sigma);
		ldf.push_back(phi);
		multiplicity.push_back(multicount/Nw);	
	}
	
	//write out the result
	FILE *result,*position;
	result = fopen("basicOutput","w");
	fprintf(result,"t	ldf	mean	var	multiplicity\n");
	for(i=0;i<Ntint;i++){
		fprintf(result,"%lf	%lf	%lf	%lf	%lf\n",(i+1)*tint,ldf[i],mean[i],var[i],multiplicity[i]);
	}
	position = fopen("position","w");
	for(i=0;i<Nw;i++){
		fprintf(position,"%lf\n",walker[i].x);
	}
	return 0;
}
