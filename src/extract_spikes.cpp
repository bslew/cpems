/*!
 * \file remove_spikes.cpp
 *
 *  Created on: Sep 8, 2010
 *      Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <tclap/CmdLine.h>
#include "Mscs-function.h"

#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#endif
//using namespace MATPACK; // use all classes and functions in the MATPACK namespace
using namespace TCLAP;

void parseOptions(int argc, char** argv);

bool _testName, _startFromAsIdx, _loop, _removeInLoop, _endAtSet, _shiftMode, _maskMode;
long _window, _shiftSize, _direction;
double _nsigma, _startFrom,_width,_endAt;
string _outfile, _infile, _removeRangeAbove, _removeRange;

//void printUsage() {
//	printf("USAGE: remove_spikes XYfile windowSize shiftSize nsigma startFrom [outfile]\n\n "
//			"where:\n\n "
//			"windowSize is the number of points in the window for calculating stdev\n"
//			"shiftSize - is the number of points by which to shift the window\n"
//			"nsigma - is the number of standard deviations defining threshold above which the points will be removed.\n"
//			"startFrom - the point number the start with (default: 0)\n"
//			"direction [12/21] - the direction of window shift 12 - increasing, 21 - decreasing (default: 0)\n"
//			"outfile - the name of the output file\n\n"
//			"when used with no arguments it will run an example and show this message.\n\n"
//			"The name of the ofputfile is given automatically based on the input file name and the parameters given");	
//}

int main(int argc, char** argv) {
	char tmpch[1000];
	mscsFunction f("f"),c("copy");
	mscsFunction spikes("spikes");
	mscsFunction extra("extra",Zero);
	long n,m,d;
	double startFrom,nsigma,endAt;

	
	parseOptions(argc, argv);

	n=_window;
	m=_shiftSize;
	nsigma=_nsigma;
	d=_direction;
	startFrom=_startFrom;
	
	if (_testName) {
				
		c=f.mkSin(0,1,0.0001,0.01);
		for (long i=0;i<f.pointsCount();i++) { f.f(i)+=10; 	i+=500; }
		f.save("input_sin_signal");
		n=500;
		m=250;
		nsigma=3;
		f.extractSpikes(n,m,nsigma);
		f.save("input_sin_signal-no-spikes");
		exit(0);
	}

	f.load(_infile);
	if (_shiftMode) {
		if (_endAtSet==false)
			spikes=f.extractSpikes(n,m,nsigma,startFrom,d,_startFromAsIdx);
		else
			spikes=f.extractSpikesRange(n,m,nsigma,startFrom,_endAt,d,_startFromAsIdx);
		c=spikes;
		printf("-- removed %li points\n",spikes.pointsCount());
		if (_loop) {
			while (c.pointsCount()>0) {
				if (_endAtSet==false)
					c=f.extractSpikes(n,m,nsigma,startFrom,d,_startFromAsIdx);
				else
					c=f.extractSpikesRange(n,m,nsigma,startFrom,_endAt,d,_startFromAsIdx);
				printf("-- removed %li points\n",spikes.pointsCount());
				spikes.concatenate(c);

				if (_width>0 and _removeInLoop) {
					mscsFunction extra("extra",Zero);
					long i;
					for (long spike = 0; spike < c.pointsCount(); spike++) {
						i=0;
						while (i<f.pointsCount()) {
							if (fabs(f.getX(i)-c.getX(spike)) < _width) { extra.newPoint(f[i]); f.deletePoint(i); }
							else i++;
						}
					}
					spikes.concatenate(extra);
				}
			}
		}
		
		if (_width>0) {
			printf("removing data around the selected spikes within range: %lE\n",_width);
			long i;
			for (long spike = 0; spike < spikes.pointsCount(); spike++) {
				i=0;
				while (i<f.pointsCount()) {
					if (fabs(f.getX(i)-spikes.getX(spike)) < _width) { extra.newPoint(f[i]); f.deletePoint(i); }
					else i++;
				}
			}
			printf("done\n");
			printf("concatenating\n");
			spikes.concatenate(extra);
			printf("done\n");
//			f.remove(extra);
		}
	}
	
	if (_maskMode) {
		
		if (_removeRangeAbove!="") {
			matrix<double> m=cpeds_matrix_load(_removeRangeAbove);
			long i;
			double x;
			extra.clearFunction();
			for (long range = 0; range < m.RowNo(); range++) {
				printf("processing range: %li of %li\n",range+1,m.RowNo());
				i=0;
				while (i<f.pointsCount()) {
					x=f.getX(i);
					if (x>m(range,0) and x<m(range,1) and f.f(i)>m(range,2)) { extra.newPoint(f[i]); f.deletePoint(i); }
					else i++;
				}
			}
			spikes.concatenate(extra);
		}

		if (_removeRange!="") {
			matrix<double> m=cpeds_matrix_load(_removeRange);
			long i;
			double x;
			extra.clearFunction();
			for (long range = 0; range < m.RowNo(); range++) {
				printf("processing range: %li of %li\n",range+1,m.RowNo());
				i=0;
				while (i<f.pointsCount()) {
					x=f.getX(i);
					if (x>m(range,0) and x<m(range,1) ) { extra.newPoint(f[i]); f.deletePoint(i); }
					else i++;
				}
			}
			spikes.concatenate(extra);
		}
	}
		
	if (_outfile=="") sprintf(tmpch,"%s-nospikes-window_%li-shift_%li-nsigma_%li",_infile.c_str(),n,m,nsigma);
	else strcpy(tmpch,_outfile.c_str());
	printf("saving signal\n");
	f.save(tmpch);
	
	if (_outfile=="") sprintf(tmpch,"%s-spikes-window_%li-shift_%li-nsigma_%li",_infile.c_str(),n,m,nsigma);
	else sprintf(tmpch,"%s.spikes",_outfile.c_str());
	printf("saving spikes\n");
	spikes.save(tmpch);		
	
	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	string::size_type j;
	
	try {
		
		CmdLine cmd("extract_spikes: spikes extractor from 2 column data files based on provided threshold in terms of number of standard deviations in a"
				"given window.\n", ' ', "" );
		
		//
		// Define arguments
		//
		
		
		UnlabeledValueArg<string> input_file("input_file","file XY function file",true,"nofile","string");	cmd.add( input_file );
		
		ValueArg<string> output("o","out","outfile ",false,"out","string"); cmd.add(output);
		ValueArg<string> removeRangeAbove("","removeRangeAbove","name of the file containing 3 columns from, to, level:\n where\n"
				"from - beginning of the window in which the data are to be removed\n"
				"to - end of the window\n"
				"level - level above which the data points will be removed\n",false,"","string"); cmd.add(removeRangeAbove);
		
		ValueArg<string> removeRange("","removeRange","name of the file containing 2 columns from, to:\n where\n"
				"from - beginning of the window in which the data are to be removed\n"
				"to - end of the window\n",false,"","string"); cmd.add(removeRange);

				
		ValueArg<long> window("n","window","Number of points in the window (default: 1000)",false,1000,"long"); cmd.add(window);
		ValueArg<long> shift("s","shift","Number of points by which to shift the window (default: 10)",false,10,"long"); cmd.add(shift);
		SwitchArg shiftMode("","shiftMode", "enable removing spikes in a mode in which a moving window passes through the signal with parameters defined"
				"by n,s,d,ns f and t.", false);	cmd.add(shiftMode);
		SwitchArg maskMode("","maskMode", "enable removing spikes in a mode in which ranges of data are defined by removeRange parameter", false);	cmd.add(maskMode);
		
		vector<long> allowedDirs;
		allowedDirs.push_back(12);
		allowedDirs.push_back(21);
		ValuesConstraint<long>* allowedDirsConstraint = new ValuesConstraint<long>(allowedDirs);
		ValueArg<long> direction("d","direction","shift direction (default:12)",false,12,allowedDirsConstraint); cmd.add(direction);
		ValueArg<double> nsigma("","ns","number of standard deviations treated as thresholf for points removal (default 3)",false,3,"double"); cmd.add(nsigma);
		ValueArg<double> startFrom("f","startFrom","In the asIds=true mode the flag indicates the point index to start from (from the left for direction 12 and from the right for direction 21)."
				"In the asIdx=false mode the flag indicates the value of the argument closest to the given number. The window will start/end on this number and will be  shifted to the right/left depending on the value of the direction flag 12/21. "
				"If not set, the full range will be considered"
				"(default 0)",false,0,"double"); cmd.add(startFrom);
		ValueArg<double> endAt("t","endAt","In the asIds=true mode the flag indicates the point index to end at (from the left for direction 12 and from the right for direction 21)."
				"In the asIdx=false mode the flag indicates the value of the argument closest to the given number. The window will start/end on this number and will be  shifted to the right/left depending on the value of the direction flag 12/21. "
						"(default 0)",false,0,"double"); cmd.add(endAt);
		SwitchArg startFromAsIdx("","asIdx", "whether or not to treat the startFrom parameter as the point index or the value of the X argument of the function", false);	cmd.add(startFromAsIdx);
				
		SwitchArg test("","test", "starts a demo test of the program.", false);	cmd.add(test);
				
		SwitchArg loop("","loop", "loops the spikes removal untill no points are removed for the current window and shift parameters."
				"In case of the looped run, the output file with spikes will contain the removed points only from the last pass of the loop", false);	cmd.add(loop);

		SwitchArg removeInLoop("","removeInLoop", "Remove window of width w also at every pass over the loop (relevant only for loop calculations).", false);	cmd.add(removeInLoop);

				
		ValueArg<double> width("w","width","The half width of the window centered on the detected spike position within which the data are to be considered as contaminated and treated as spikes. "
				"For 0 (default 0) value only the detected spike will be in the output file.",false,0,"double"); cmd.add(width);
				
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_infile = input_file.getValue(); /*if ( _input_file == "nofile" || _input_file == "Ylm" ) { _not_from_file = true; } else { _not_from_file = false; }*/
		_outfile = output.getValue();
		_removeRangeAbove= removeRangeAbove.getValue();
		_removeRange= removeRange.getValue();
		_window= window.getValue();
		_shiftSize = shift.getValue();
		_direction = direction.getValue();
		// _loffset = loffset.getValue();
		_nsigma = nsigma.getValue();
		_startFrom = startFrom.getValue();
		_endAt = endAt.getValue(); if (endAt.isSet()) _endAtSet=true; else _endAtSet=false; 
		_startFromAsIdx = startFromAsIdx.getValue();
		_testName = test.getValue();
		_width=width.getValue();
		_shiftMode=shiftMode.getValue();
		_maskMode=maskMode.getValue();
		
		_loop= loop.getValue();
		_removeInLoop= removeInLoop.getValue();
				
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



