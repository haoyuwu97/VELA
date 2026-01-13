#include "correlator.h"
#include <math.h>
#include <algorithm>
#include <vector>

#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

void Correlator::initialize_Gw() {

	for (unsigned int i=0;i<length;++i) {
		wt[i] = 0;
		Gw[i] = 0;
		Gw_storage[i] = 0;
		Gw_loss[i] = 0;
	}

}

void Correlator::calculate_Gw(const double dtime, const int time, const int Gtcut, const double tail_param) {
	
	double X1=0.0,X2=0.0,costk_1=0.0,sintk_1=0.0,coswdtk=0.0,sinwdtk=0.0,dtk,w2;
	double dt=dtime;
	int tt=time;
	double* wtt=new double [tt];
	double* Gw_1=new double [tt];
	double* Gw_2=new double [tt];
	unsigned int im=0;
	int max_index = (Gtcut > 0 && Gtcut < static_cast<int>(npcorr)) ? Gtcut : static_cast<int>(npcorr);
	std::vector<double> Gt_work(npcorr, 0.0);

	for (unsigned int i=0;i<tt;++i) {
		wtt[i] = 0.0;
		Gw_1[i] = 0.0;
		Gw_2[i] = 0.0;
	}

	if (tt <= 1 || npcorr == 0) {
		delete [] wtt;
		delete [] Gw_1;
		delete [] Gw_2;
		return;
	}

	for (int i=0;i<max_index;++i) {
		Gt_work[i] = Gt[i];
	}

	int window = std::max(5, max_index / 10);
	int consecutive = std::max(3, window / 4);
	double mean = 0.0;
	double var = 0.0;
	int tail_start = std::max(0, max_index - window);
	for (int i=tail_start;i<max_index;++i) {
		mean += Gt_work[i];
	}
	mean /= static_cast<double>(max_index - tail_start);
	for (int i=tail_start;i<max_index;++i) {
		double diff = Gt_work[i] - mean;
		var += diff * diff;
	}
	var /= static_cast<double>(std::max(1, max_index - tail_start - 1));
	double sigma = sqrt(var);
	double k = tail_param > 0.0 ? tail_param : 2.0;
	int run = 0;
	int cutoff = max_index;
	for (int i=max_index-1;i>=0;--i) {
		if (fabs(Gt_work[i]) <= k * sigma) {
			run++;
			if (run >= consecutive) {
				cutoff = i;
				break;
			}
		} else {
			run = 0;
		}
	}
	if (cutoff < max_index) {
		max_index = std::max(cutoff, 2);
	}

	for (unsigned int i=1;i<tt;++i) {
		wtt[i]=1/(i*dt);
		w2=wtt[i]*wtt[i];
		for (int j=1;j<max_index;++j){
			
			costk_1 = cos(wtt[i]*t[j-1]*dt);
			sintk_1 = sin(wtt[i]*t[j-1]*dt);
			dtk = (t[j]-t[j-1])*dt;
			coswdtk = cos(wtt[i]*dtk);
			sinwdtk = sin(wtt[i]*dtk);
			
			X1 = Gt_work[j]*( (coswdtk-1)/(w2*dtk) + sinwdtk/wtt[i] ) - Gt_work[j-1]*( (coswdtk-1)/(w2*dtk) );
			X2 = Gt_work[j]*( coswdtk/wtt[i] - sinwdtk/(w2*dtk) ) - Gt_work[j-1]*( 1.0/wtt[i] - sinwdtk/(w2*dtk) );
			
			Gw_1[i] += sintk_1*X1 - costk_1*X2;
			Gw_2[i] += sintk_1*X2 + costk_1*X1;
		}
		Gw_1[i] = Gw_1[i]*wtt[i];
		Gw_2[i] = Gw_2[i]*wtt[i];
		//std::cout << wtt[i] << " " << Gw_1[i] << " " << Gw_2[i] << std::endl;
	}

	for(int i=0;i<p;++i){
		wt[im]=i;
		Gw_storage[im]=Gw_1[im];
		Gw_loss[im]=Gw_2[im];
		Gw[im] = sqrt( Gw_storage[im]*Gw_storage[im] + Gw_loss[im]*Gw_loss[im] );
		++im;
	}
	int ss=im;
	for (int k=1;k<kmax;++k) {
		for (int i=dmin;i<p;++i) {
			wt[im]=i * pow((double)m, k);
			++im;
		}
	}

	for (int k=1;k<kmax;++k) {
		for (int i=dmin;i<p;++i) {
			int lgth=0;
			if(ss+1>=npcorr) break;
			int start = static_cast<int>(wt[ss-1]);
			int end = static_cast<int>(wt[ss+1]);
			start = std::max(0, start);
			end = std::min(end, tt - 1);
			if (start > end) {
				continue;
			}
			for(int s=start;s<=end;++s){
				Gw_storage[ss]+=Gw_1[s];
				Gw_loss[ss]+=Gw_2[s];
				lgth++;
			}
			Gw_storage[ss]=Gw_storage[ss]/lgth;
			Gw_loss[ss]=Gw_loss[ss]/lgth;
			Gw[ss] = sqrt( Gw_storage[ss]*Gw_storage[ss] + Gw_loss[ss]*Gw_loss[ss] );
			ss++;
		}
	}

	delete [] wtt;
	delete [] Gw_1;
	delete [] Gw_2;

}

