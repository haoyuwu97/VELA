#include "correlator.h"
#include <math.h>
#include <cmath>

#include <iostream>
#include <fstream>
#include <stdlib.h>

int main(int argc, char *argv[]) {

	if (argc < 6) {
		std::cout << "Usage:" << std::endl 
			<< "Program INPUT_FILE OUTPUT_FILE_G(T) OUTPUT_FILE_G(W) DT(YOUR_SIMALLESTDATA=DT*TIME_UNIT) MODE [G(T)CUT] [TAIL_PARAM]" << std::endl
			<< "Note: You need to multiply the result of the modulus by V/kBT" << std::endl
			<< "      MODE==0:G(t) used stress in xy,xz,yz" << std::endl
			<< "      MODE==1:G(t) used stress in xx,yy,zz,xy,xz,yz; and your input file should put xx,yy,zz in first three columns" << std::endl
			<< "      G(t)cut sets a maximum time index for tail handling; defaults to the time maximum" << std::endl
			<< "      TAIL_PARAM is the k-sigma threshold for the statistical confidence cutoff" << std::endl;
		return 1;
	}

	Correlator c;
	c.setsize(32, 16, 2);    

	std::ifstream fin;
	double valab,valac,valbc,valaa,valbb,valcc;
	
	fin.open(argv[1]);
	std::ofstream fout;
	fout.open(argv[2]);
	std::ofstream fout2;
	fout2.open(argv[3]);
	double dt=atof(argv[4]);
	int mode=atoi(argv[5]);
	int Gtcut=-1;
	double tail_param=2.0;
	if (argc >= 7){
		Gtcut=atoi(argv[6]);
		if(Gtcut<0){
			std::cout << "ERROR G(T)CUT!" << std::endl;
			return 1;
		}
	}
	if (argc >= 8){
		tail_param=atof(argv[7]);
	}
	if(!fin.is_open() || fin.bad()) {
		std::cout << "ERROR INPUT FILE!" << std::endl;
		return 1;
	}
	if(!fout.is_open() || fout.bad() || !fout2.is_open() || fout2.bad()) {
		std::cout << "ERROR OUTPUT FILE!" << std::endl;
		return 1;
	}
	if (dt <= 0.0) {
		std::cout << "ERROR DT!" << std::endl;
		return 1;
	}
	if(mode!=0 && mode!=1) {
		std::cout << "ERROR MODE SET!" << std::endl;
		return 1;
	}
	
	int tt;
	long line_number = 0;
	if(mode==0){
		c.initialize_mode0();
		tt=0;
		while(fin >> valab >> valac >> valbc) {
			++line_number;
			if (!std::isfinite(valab) || !std::isfinite(valac) || !std::isfinite(valbc)) {
				std::cout << "ERROR NON-FINITE INPUT AT LINE " << line_number << "!" << std::endl;
				return 1;
			}
			c.add_mode0 (valab,valac,valbc);
			tt++;
		}
		if (!fin.eof()) {
			std::cout << "ERROR INVALID INPUT FORMAT AT LINE " << (line_number + 1) << "!" << std::endl;
			return 1;
		}
		fin.close();
		if (tt == 0) {
			std::cout << "ERROR EMPTY INPUT!" << std::endl;
			return 1;
		}
		c.evaluate_mode0();
		fout << "#dt=" << dt << " mode=" << mode << " Gtcut=" << Gtcut
			<< " tail_param=" << tail_param << std::endl;
		fout << "#Time " << "G_t " << std::endl;
		for (unsigned int i=0;i<c.npcorr;++i){
			c.Gt[i] = (c.f1_mode0[i] + c.f2_mode0[i] + c.f3_mode0[i]) / 3.0;
			fout << c.t[i]*dt << " " << c.Gt[i] << std::endl;
		}
		fout.close();
	}else if(mode==1){
		c.initialize_mode1();
		tt=0;
		while(fin >> valaa >> valbb >> valcc >> valab >> valac >> valbc) {
			++line_number;
			if (!std::isfinite(valaa) || !std::isfinite(valbb) || !std::isfinite(valcc) ||
				!std::isfinite(valab) || !std::isfinite(valac) || !std::isfinite(valbc)) {
				std::cout << "ERROR NON-FINITE INPUT AT LINE " << line_number << "!" << std::endl;
				return 1;
			}
			//c.add_mode1 (valaa-valbb,valaa-valcc,valbb-valcc,valab,valac,valbc);
			c.add_mode1 (valaa,valbb,valcc,valab,valac,valbc);
			tt++;
		}
		if (!fin.eof()) {
			std::cout << "ERROR INVALID INPUT FORMAT AT LINE " << (line_number + 1) << "!" << std::endl;
			return 1;
		}
		fin.close();
		if (tt == 0) {
			std::cout << "ERROR EMPTY INPUT!" << std::endl;
			return 1;
		}
		c.evaluate_mode1();
		fout << "#dt=" << dt << " mode=" << mode << " Gtcut=" << Gtcut
			<< " tail_param=" << tail_param << std::endl;
		fout << "#Time " << "G_t " << std::endl;
		for (unsigned int i=0;i<c.npcorr;++i){
			c.Gt[i] = (c.f1_mode1[i] + c.f2_mode1[i] + c.f3_mode1[i]) / 30.0 + (c.f4_mode1[i] + c.f5_mode1[i] + c.f6_mode1[i]) / 5.0;
			fout << c.t[i]*dt << " " << c.Gt[i] << std::endl;
		}
		fout.close();
	}
	
	c.initialize_Gw();
	if(Gtcut==-1) Gtcut=c.npcorr;
	c.calculate_Gw(dt,tt,Gtcut,tail_param);
	fout2 << "#dt=" << dt << " mode=" << mode << " Gtcut=" << Gtcut
		<< " tail_param=" << tail_param << std::endl;
	fout2 << "#Frequency(1/time_unit) " << "Gw_storage " << "Gw_loss " << "Gw " << std::endl;
	for (unsigned int i=c.npcorr-1;i>0;--i){
		fout2 << 1.0/(c.wt[i]*dt) << " " << c.Gw_storage[i] << " " << c.Gw_loss[i] << " " << c.Gw[i] << std::endl;
	}
	fout2.close();

	return 0;
}
